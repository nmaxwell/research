

class ccFloat:
    
    value_=0.0
    
    def __init__(self, rhs=0.0):
        self.value=rhs
    
    
    def __add__(self, rhs):
        return self.value.__add__(rhs)
    
    def __add__(self, rhs):
        return self.value.__add__(rhs)
    
    
    
    
    def __lt__(self, rhs):
        return self.value < rhs
    
    def __le__(self, rhs):
        return self.value <+ rhs
    
    def __eq__(self, rhs):
        return self.value == rhs
    
    def __ne__(self, rhs):
        return self.value != rhs
    
    def __gt__(self, rhs):
        return self.value > rhs
    
    def __lt__(self, rhs):
        return self.value >= rhs
    
    
    
    
    
    
    
    def __str__(self ):
        return str(self.value)
    
    
    










