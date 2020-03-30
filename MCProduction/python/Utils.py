"""_Utils_

Module containing some utility tools

"""

def useSystemRandomSeeds(process):
    """_useSystemRandomSeeds_
    
    Initiate RandomNumberServiceHelper seed using random.SystemRandom
    to have different seeds every time one runs cmsRun.

    """
    from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
    randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService) 
    randSvc.populate()

    return process

def registerAndParseArgument(name):
    """_registerAndParseArgument_
    
    Register and parse a new command line argument. Return VarParsing object.
    """
    from FWCore.ParameterSet.VarParsing import VarParsing
    options = VarParsing ('analysis')
    options.register(name, 0, VarParsing.multiplicity.singleton,VarParsing.varType.int,name)
    options.parseArguments()
    return options
