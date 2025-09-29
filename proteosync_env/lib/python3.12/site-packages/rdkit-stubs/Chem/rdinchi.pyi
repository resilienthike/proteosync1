from __future__ import annotations
__all__: list[str] = ['GetInchiVersion', 'InchiToInchiKey', 'InchiToMol', 'MolBlockToInchi', 'MolToInchi', 'MolToInchiKey']
def GetInchiVersion() -> str:
    """
        returns the version of the InChI software being used
    
        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetInchiVersion()
    """
def InchiToInchiKey(inchi: str) -> str:
    """
        return the InChI key for an InChI string
    
        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > InchiToInchiKey(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def InchiToMol(inchi: str, sanitize: bool = True, removeHs: bool = True) -> tuple:
    """
        return a ROMol for a InChI string
          Returns:
            a tuple with:
              - the molecule
              - the return code from the InChI conversion
              - a string with any messages from the InChI conversion
              - a string with any log messages from the InChI conversion
        
    
        C++ signature :
            boost::python::tuple InchiToMol(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True]])
    """
def MolBlockToInchi(molblock: str, options: str = '') -> tuple:
    """
        return the InChI for a ROMol molecule.
        
          Arguments:
            - molblock: the mol block to use.
            - options: the InChI generation options.
              Options should be prefixed with either a - or a /
              Available options are explained in the InChI technical FAQ:
              https://www.inchi-trust.org/technical-faq/#15.14
              and the User Guide:
              https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
          Returns:
            a tuple with:
              - the InChI
              - the return code from the InChI conversion
              - a string with any messages from the InChI conversion
              - a string with any log messages from the InChI conversion
              - a string with the InChI AuxInfo
        
    
        C++ signature :
            boost::python::tuple MolBlockToInchi(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=''])
    """
def MolToInchi(mol: Mol, options: str = '') -> tuple:
    """
        return the InChI for a ROMol molecule.
        
          Arguments:
            - mol: the molecule to use.
            - options: the InChI generation options.
              Options should be prefixed with either a - or a /
              Available options are explained in the InChI technical FAQ:
              https://www.inchi-trust.org/technical-faq/#15.14
              and the User Guide:
              https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
          Returns:
            a tuple with:
              - the InChI
              - the return code from the InChI conversion
              - a string with any messages from the InChI conversion
              - a string with any log messages from the InChI conversion
              - a string with the InChI AuxInfo
        
    
        C++ signature :
            boost::python::tuple MolToInchi(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=''])
    """
def MolToInchiKey(mol: Mol, options: str = '') -> str:
    """
        return the InChI key for a ROMol molecule.
        
          Arguments:
            - mol: the molecule to use.
            - options: the InChI generation options.
              Options should be prefixed with either a - or a /
              Available options are explained in the InChI technical FAQ:
              https://www.inchi-trust.org/technical-faq/#15.14
              and the User Guide available from:
              https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
          Returns: the InChI key
        
    
        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToInchiKey(RDKit::ROMol [,char const*=''])
    """
