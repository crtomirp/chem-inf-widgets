#!/usr/bin/env python
import argparse
import numpy as np
from rdkit import Chem
import torch

class SmilesEnumerator:
    """
    SMILES Enumerator, vectorizer and devectorizer.
    
    Arguments:
        charset: A string containing the characters for vectorization.
        pad: Length to pad the vectorization.
        leftpad: If True, pads on the left side.
        isomericSmiles: Include stereochemistry information.
        enum: Enumerate (randomize) the SMILES during transform.
        canonical: Use canonical SMILES (overrides enumeration).
    """
    def __init__(self, charset='@C)(=cOn1S2/H[N]\\', pad=120, leftpad=True, 
                 isomericSmiles=True, enum=True, canonical=False):
        self._charset = None
        self.charset = charset
        self.pad = pad
        self.leftpad = leftpad
        self.isomericSmiles = isomericSmiles
        self.enumerate = enum
        self.canonical = canonical

    @property
    def charset(self):
        return self._charset
        
    @charset.setter
    def charset(self, charset):
        self._charset = charset
        self._charlen = len(charset)
        self._char_to_int = {c: i for i, c in enumerate(charset)}
        self._int_to_char = {i: c for i, c in enumerate(charset)}
        
    def fit(self, smiles, extra_chars=[], extra_pad=5):
        """
        Extracts the charset from the SMILES dataset and sets the pad length.
        
        Arguments:
            smiles: Array-like of SMILES strings.
            extra_chars: List of extra characters to add to the charset.
            extra_pad: Extra padding to add to the length.
        """
        charset = set("".join(smiles))
        self.charset = "".join(charset.union(set(extra_chars)))
        self.pad = max(len(s) for s in smiles) + extra_pad
        
    def randomize_smiles(self, smiles):
        """
        Performs a randomization of a SMILES string (if the molecule is valid).
        """
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            raise ValueError("Invalid SMILES string: " + smiles)
        atom_indices = list(range(m.GetNumAtoms()))
        np.random.shuffle(atom_indices)
        nm = Chem.RenumberAtoms(m, atom_indices)
        return Chem.MolToSmiles(nm, canonical=self.canonical, 
                                isomericSmiles=self.isomericSmiles)

    def transform(self, smiles):
        """
        Enumerates (randomizes) and vectorizes a single SMILES string into a one-hot tensor.
        
        Arguments:
            smiles: A single SMILES string.
            
        Returns:
            one_hot: A torch.Tensor of shape (pad, _charlen).
        """
        one_hot = torch.zeros((self.pad, self._charlen), dtype=torch.int8)
        s = smiles
        if self.enumerate:
            s = self.randomize_smiles(s)
        if self.leftpad:
            l = len(s)
            diff = self.pad - l
            for j, c in enumerate(s):
                one_hot[j + diff, self._char_to_int[c]] = 1
        else:
            for j, c in enumerate(s):
                one_hot[j, self._char_to_int[c]] = 1
        return one_hot

    def reverse_transform(self, vect):
        """
        Converts a one-hot encoded tensor back to a SMILES string.
        
        Arguments:
            vect: A 2D torch.Tensor of shape (pad, _charlen).
            
        Returns:
            A SMILES string.
        """
        mask = vect.sum(dim=1) == 1
        indices = torch.argmax(vect[mask], dim=1)
        smile = "".join(self._int_to_char[int(i)] for i in indices)
        return smile

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate random SMILES codes for caffeine")
    parser.add_argument("--num", type=int, default=100, help="Number of random SMILES to generate")
    args = parser.parse_args()
    
    # Caffeine canonical SMILES string (the starting structure)
    caffeine_smiles = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
    
    # Instantiate the enumerator with randomization enabled.
    enumerator = SmilesEnumerator(canonical=False, enum=True)
    
    # Generate and print the requested number of random SMILES codes.
    for i in range(args.num):
        try:
            randomized_smile = enumerator.randomize_smiles(caffeine_smiles)
            print(randomized_smile)
        except Exception as e:
            print(f"Error generating SMILES: {e}")
