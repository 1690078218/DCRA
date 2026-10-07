from .main import main as _pipeline

def dcra(argv = None):
    if argv == None:
        _pipeline('-h')
        return 
    else: 
        return _pipeline(argv)

