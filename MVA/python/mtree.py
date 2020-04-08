from ROOT import TTree
from array import array
import re

class MTree:
    def __init__(self, name, title):
        self.tree = TTree(name,title)
        # map between TTree types and array types
        # https://root.cern.ch/root/html524/TTree.html
        # https://docs.python.org/2/library/array.html
        self.known_types = {
            'Int_t':{'array':'i','root':'I'}, 
            'UInt_t':{'array':'I','root':'i'}, 
            'Float_t':{'array':'f','root':'F'},
            'ULong64_t':{'array':'L','root':'l'}
        }
        self.defaults = dict()
        self.variables = dict()

    def addBranch(self, branch_name, branch_type, default_value, title = None):
        """Register new branch.

        The new branch can be of one of the following types: Int_t,
        UInt_t, Float_t and ULong64_t. The default value is used to
        later to reset branch values to their defaults. Title is
        optional, but highly advisable since it helps to make the
        ntuple self documenting.
        """
        if branch_type not in self.known_types:
            raise Exception("Uknown type %s" % branch_type)
        array_type = self.known_types[branch_type]['array']
        root_type = self.known_types[branch_type]['root']
        self.defaults[branch_name] = default_value
        self.variables[branch_name] = array(array_type,[default_value])
        # setattr(self, branch_name, self.variables[branch_name][0])
        self.tree.Branch(branch_name, self.variables[branch_name], "%s/%s"%(branch_name,root_type))
        if title:
            self.tree.GetBranch(branch_name).SetTitle(title)
        
    def __setitem__(self, branch_name, value):
        self.variables[branch_name][0] = value

    def __getitem__(self, branch_name):
        return self, branch_name

    def reset(self, branch_list=None, regexp=None):
        """Initialize variables to their default values.

        In order to avoid storing incorrect information it is
        advisable to reset branches to their default values. One may
        restrict for branches to reset using branch_list or regular
        expressions.
        """
        for branch_name,value in self.defaults.items():
            if branch_list == None or branch_name in branch_list:
                if regexp == None or re.search(regexp,branch_name):
                    self.variables[branch_name][0] = value

    def fill(self):
        """Store current branch values"""
        self.tree.Fill()

def main():
    t = MTree("test","")
    t.addBranch("pt","Float_t",0)
    t.addBranch("sim_pt","Float_t",0)
    t.addBranch("gen_evt_pt","Float_t",0)
    t.addBranch("evt_pt","Float_t",0)

    t["pt"] = 3.1
    t["sim_pt"] = 3.2
    t.fill()

    t.reset()
    t["pt"] = 14.8
    t["sim_pt"] = 14.3
    t.fill()

    t.reset(regexp="^(?!pt).")
    t["evt_pt"] = 2.7

    t.fill()

    t.reset(regexp="pt")
    t.fill()

    t.reset()
    t['gen_evt_pt'] = 13.1
    t['evt_pt'] = 14.2

    t.fill()

    t.reset(regexp="^(?!evt_).")
    t.fill()

    t.tree.Print()
    
    t.tree.Scan()
    

if __name__ == "__main__":
    # help(MTree)
    main()
