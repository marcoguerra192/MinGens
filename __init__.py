import sys,os

MinGensDir = os.path.abspath(os.path.dirname(__file__))

if MinGensDir not in sys.path:
    
    sys.path.append(MinGensDir)