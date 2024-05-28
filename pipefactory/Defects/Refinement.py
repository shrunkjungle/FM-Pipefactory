class AxialRefinement:
    def __init__(self, distance, new_element_size, func):
        self.x = distance
        self.dx = new_element_size
        self.func = func

    def __call__(self, length):
        return self.func(length - self.x)
    
class Ramp:
    def __init__(self, w0, w1):
        """
        Defines a linear ramp class with function decreasing from 1 to 0 from x =`w0` to x =`w1` respectively
        """
        self.w0 = w0
        self.w1 = w1
    
    def __call__(self, inp):
        x = abs(inp)
        if x <= self.w0:
            return 1.0
        elif x > self.w0 and x<self.w1:
            return (self.w1 - x)/(self.w1-self.w0)
        else:
            return 0.0