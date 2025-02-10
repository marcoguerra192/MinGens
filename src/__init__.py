import sys,os

srcDir = os.path.abspath(os.path.dirname(__file__))

if srcDir not in sys.path:
	sys.path.append(srcDir)
