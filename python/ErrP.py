# ErrP class definition


class ErrP:
    
    def __init__(self,name='MyModel',params=None,parunc=None,covmatrix=None,func_model=None):
        """Constructor. Optional parameter name."""
        
        #print('Creating object')
        
        # informations of model
        self.name=name
        
        self.params = params        # model parameters
        self.parunc = parunc        # model parameters uncertainties
        self.covmarix = covmatrix   # covariance matrix
        self.funcmodel = func_model # model function
        
        
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
    
    def setModelName(self,name):
        """Set model name"""
        
        self.name=name
        
    def useCovMatrix(option=False):
        """Set if user wants to use covariance matrix"""
        self.usecov = option
        
    def setModel(self,func_model):
        """Set model function"""
        self.funcmodel = func_model
        
    def evalModel(self,x):
        """Evaluates the model passed by the user with given parameters for x"""
        model = self.funcmodel
        par = self.params
        if self.parstatus == True: #not enough to solve the problem
            return model(x,par)
        else:
            return model(x)
    