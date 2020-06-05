# ErrP class definition

import numpy as np #functions need numpy module

class ErrP:

    def __init__(self,name='MyModel',
                 params=None,
                 parunc=None,
                 covmatrix=None,
                 func_model=None):

        """Constructor. Optional parameter name."""

        #print('Creating object')

        # informations of model
        self.name=name

        self.params = params        # model parameters
        self.parunc = parunc        # model parameters uncertainties
        self.covmatrix = covmatrix   # covariance matrix
        self.model = func_model # model function


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
        """Set model name"""
        self.name=name

    def UseCovMatrix(self,option=False):
        """Set if user wants to use covariance matrix"""
        self.usecov = option

    def SetModel(self,func_model):
        """Set model function"""
        self.model = func_model

    def SetParameters(self,user_par):
        self.params = user_par
        self.parstatus = True

    def SetParErrors(self,user_errors):
        """Set parameters uncertantines"""
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


    def CalcModel(self,x):
        """Evaluates the model passed by the user with given parameters for x"""
        model = self.model
        par = self.params
        #function expect parameters
        #if self.parstatus == True: #not enough to solve the problem
        return model(x,par)
        #else:
         #   return model(x)


    def CalcGrad(self,x): #calculate gradient
        """ Evaluates the Gradient of Model """

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

        params = self.params
        sigma = self.unc

        parmc = np.random.normal(params,sigma)

        return parmc

    def GenParMCCov(self,x):

        params = self.params
        covmatrix = self.covmatrix

        parmc = np.random.multivariate_normal(params,covmatrix)

        return parmc



    #def evalMCError(self,x) # evaluate error prop with MC
