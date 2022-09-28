import sys
import numpy as np


''''''
def PAMParser(pam_path: str) -> list:
    pam = []

    with open(pam_path, 'r') as pam:
        line = pam.readline()
