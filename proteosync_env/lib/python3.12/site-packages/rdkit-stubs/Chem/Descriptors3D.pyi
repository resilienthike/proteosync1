"""
 Descriptors derived from a molecule's 3D structure

"""
from __future__ import annotations
from rdkit.Chem.Descriptors import _isCallable
from rdkit.Chem import rdMolDescriptors
__all__: list[str] = ['CalcMolDescriptors3D', 'descList', 'rdMolDescriptors']
def CalcMolDescriptors3D(mol, confId = None):
    """
    
        Compute all 3D descriptors of a molecule
        
        Arguments:
        - mol: the molecule to work with
        - confId: conformer ID to work with. If not specified the default (-1) is used
        
        Return:
        
        dict
            A dictionary with decriptor names as keys and the descriptor values as values
    
        raises a ValueError 
            If the molecule does not have conformers
        
    """
def _setupDescriptors(namespace):
    ...
descList: list  # value = [('PMI1', <function <lambda> at 0x7f44e882f880>), ('PMI2', <function <lambda> at 0x7f44e882ff60>), ('PMI3', <function <lambda> at 0x7f44e73e4040>), ('NPR1', <function <lambda> at 0x7f44e73e40e0>), ('NPR2', <function <lambda> at 0x7f44e73e4180>), ('RadiusOfGyration', <function <lambda> at 0x7f44e73e4220>), ('InertialShapeFactor', <function <lambda> at 0x7f44e73e42c0>), ('Eccentricity', <function <lambda> at 0x7f44e73e4360>), ('Asphericity', <function <lambda> at 0x7f44e73e4400>), ('SpherocityIndex', <function <lambda> at 0x7f44e73e44a0>), ('PBF', <function <lambda> at 0x7f44e73e4540>)]
