# ErrP class definition

import numpy as np #functions need numpy module

class ErrP:
    
    def __init__(self,name='MyModel',params=None,parunc=None,covmatrix=None,func_model=None):
        """Constructor. Optional parameter name."""
        
        #print('Creating object')
        
        # informations of model
        self.name=name
        
        self.params = params        # model parameters
        self.parunc = parunc        # model parameters uncertainties
        self.covmarix = covmatrix   # covariance matrix
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
Using Covariance: {}""".format(self.name,self.funcstatus,self.parstatus,self.uncstatus,self.covstatus,self.usecov)
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
        
        grad=[]
        par_central = self.params[:] # copy of central parameters
        par_aux = self.params[:]
        eps = 1.0e-8
        model = self.model
        
        for i in range(len(self.params)):
            par_aux[i] += eps
            delta = self.model(x,par_aux) - model(x,par_central)
            grad.append(delta/eps)
            par_aux[i] = par_central[i]
            
        return grad       
        
    #def CalcErrorGrad(self,x) #calculate error prop with derivatives
    #def genParams(self,x)   # generate normal distributed parameters
    #def evalMCError(self,x) # evaluate error prop with MC
    