from enum import Enum


class FeatureType(Enum):
    RESIDUE = "RESIDUE"
    PAIR = "PAIR"
    MSA = "MSA"
    TEMPLATES = "TEMPLATES"
    TOKEN = "TOKEN"
    TOKEN_PAIR = "TOKEN_PAIR"
    ATOM = "ATOM"
    ATOM_PAIR = "ATOM_PAIR"
