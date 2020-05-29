# ErrP class definition


class ErrP:
    
    def __init__(self,name='model'):
        """Constructor. Optional parameter name."""
        
        print('Creating object')
        self.name=name
        self.covstatus = False
        self.uncstatus = False
        self.funcstatus = False
        
    def __str__(self):
        """Conversion to string."""
        text = 'Model {}; Using Covariance: {}.'.format(self.name,self.covstatus)
        return text
    
    def setModelName(self,name):
        self.name=name