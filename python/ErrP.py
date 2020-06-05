# ErrP class definition

import numpy as np #functions need numpy module

class ErrP:

    def __init__(self,name='MyModel',
                 params=None,
                 parunc=None,
                 covmatrix=None,
                 func_model=None,
                 nsample = 1):

        """Constructor.
        Parameters (optional):
        'name': name of the model;
        'params': parameter central values;
        'parunc': parameter uncertainties;
        'covmatrix': covariance matrix (ndarray or path to file with matrix);
        'func_model': model function;
        'nsample': number of samples to be calculated when using Monte Carlo.

        All these parameters can be set using functions defined below.
        """

        #print('Creating object')

        # informations of model
        self.name=name

        self.params = params        # model parameters
        self.parunc = parunc        # model parameters uncertainties
        self.covmatrix = covmatrix   # covariance matrix
        self.model = func_model # model function
        self.nsample = nsample

        #status variables

        #status of parameters
        if params != None:
            self.parstatus = True
        else:
            self.parstatus = False

        #status of parameters uncertainties
        if parunc != None:
            self.uncstatus = True
        else:
            self.uncstatus = False

        #status of covariance matrix and if user wants to use it
        #if user informed, it is assumed that it is to be used
        if covmatrix != None:
            self.covstatus = True
            self.usecov = True
        else:
            self.covstatus = False
            self.usecov = False

        #status of model function
        if func_model != None:
            self.funcstatus = True
        else:
            self.funcstatus = False


    def __str__(self):
        """Conversion to string."""
        text = """Model {}
Function defined: {}
Parameters informed: {}
Uncertainties informed: {}
Cov matrix informed: {}
Using Covariance: {}""".format(self.name,self.funcstatus,self.parstatus,
                               self.uncstatus,self.covstatus,self.usecov)
        return text


    def SetModelName(self,name):
        """Set model name."""
        self.name=name

    def UseCovMatrix(self,option=False):
        """Set if user wants to use covariance matrix."""
        self.usecov = option

    def SetModel(self,func_model):
        """Set model function."""
        self.model = func_model

    def SetParameters(self,user_par):
        """Set parameters central values."""
        self.params = user_par
        self.parstatus = True

    def SetParErrors(self,user_errors):
        """Set parameters uncertantines."""
        self.parunc = user_errors
        self.uncstatus = True

    def SetCovMatrix(self,user_matrix):
        """Set covariance matrix.
        User matrix can be a file path or a matrix (numpy array) itself"""

        if(type(user_matrix)==str):
#            print('reading matrix from file')
            self.covmatrix = np.loadtxt(user_matrix)
        else:
#            print('matrix passed by the user')
            self.covmatrix = user_matrix

        self.covstatus = True
        self.usecov = True

    def SetNsample(self,x):
        """Set the number of samples to be generated when using MonteCarlo."""
        self.nsample = int(x)


    def CalcModel(self,x):
        """Evaluates the model passed by the user with given parameters for x."""
        model = self.model
        par = self.params
        #function expect parameters
        #if self.parstatus == True: #not enough to solve the problem
        return model(x,par)
        #else:
         #   return model(x)


    def CalcGrad(self,x): #calculate gradient
        """ Evaluates the Gradient of the model in point x.
        Derivatives are obtained with forward finite diferences.
        Returns a list with gradient components."""

        if self.funcstatus == False:
            print('Function model not informed')
            return None

        if self.parstatus == False:
            print('Central parameters not informed')
            return None

        npar = len(self.params)
        grad = [] #np.empty((npar))
        par_central = self.params[:] # copy of central parameters
        par_aux = self.params[:]
        eps = 1.0e-8
        model = self.model

        for i in range(npar):
            par_aux[i] = par_aux[i] + eps
            delta = model(x,par_aux) - model(x,par_central)
            grad.append(delta/eps)
            par_aux[i] = par_central[i]

        return grad


    def CalcErrorGrad(self,x): #calculate error prop with derivatives
        """Calculates error using derivatives (gradient),
        with or without covariance matrix.
        Returns the error.
        """

        if self.funcstatus == False:
            print('Function model not informed')
            return None
        elif self.parstatus == False:
            print('Parameters not informed')
            return None
        elif self.uncstatus == False and self.covstatus == False:
            print('Uncertainties and covariance matrix not informed')
            return None
        elif self.uncstatus == True and self.usecov == False:
            print('Error propagation without covariance')
            # set the appropriate cov matrix
            npar = len(self.params)
            unc = self.parunc
            covmatrix = np.empty((npar,npar),float)

            for i in range(npar):
                for j in range(npar):
                    if i==j:
                        covmatrix[i][j] = unc[i]**2
                    else:
                        covmatrix[i][j] = 0.0

        else:
            print('Error propagation with covariance')
            covmatrix = np.array(self.covmatrix) # converts to numpy array


        grad = np.array([self.CalcGrad(x)]) #converts list in numpy array (line matrix)
        gradT = grad.T # transpose (column matrix)

        if grad.all() == None:
            print('Error in gradient calculation. Aborting...')
            return None

        gradT_cov = covmatrix.dot(gradT)
        error2 = grad.dot(gradT_cov)[0][0]

        return np.sqrt(error2)


# generate normal distributed parameters
    def GenParMC(self,x):
        """Generate normal distributed parameters from
        central values and uncertainties (without covariance matrix).
        Returns a ndarray with parameters values.
        """

        params = self.params
        sigma = self.unc

        parmc = np.random.normal(params,sigma)

        return parmc

    def GenParMCCov(self,x):
        """Generate normal distributed parameters from
        central values and covariance matrix.
        Returns a ndarray with parameters values.
        """

        params = self.params
        covmatrix = self.covmatrix

        parmc = np.random.multivariate_normal(params,covmatrix)

        return parmc

    def CalcErrorMC(self,x):
        """Calculates Error from normal distributed random parameters generated
        from central values and uncertainties, with or without covariance
        matrix. Returns a dictionary with two keys:
        'mean': mean value of central value from Nsample samples;
        'std': standard deviavion from the Nsample generated and
               corresponds to the propagated error.
        """

        Nsample = self.nsample
        model = self.model
        vcentral = []

        if self.funcstatus == False:
            print('Function model not informed')
            return None
        elif self.parstatus == False:
            print('Parameters not informed')
            return None
        elif self.uncstatus == False and self.covstatus == False:
            print('Uncertainties and covariance matrix not informed')
            return None
        elif self.uncstatus == True and self.usecov == False:
            print('Error propagation without covariance')

            for _ in range(Nsample):
                parmc = self.GenParMC(x)
                vcentral.append(model(x,parmc))

        else:
            print('Error propagation with covariance')

            for _ in range(Nsample):
                parmc = self.GenParMCCov(x)
                vcentral.append(model(x,parmc))

        # calculating average values and standard deviations

        vcentral_array = np.array(vcentral)
        mean = np.mean(vcentral_array)
        std = np.std(vcentral_array,ddof=1)
        res = {'mean': mean, 'std': std}

        return res
