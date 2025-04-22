from __future__ import print_function
from pySecDec.integral_interface import DistevalLibrary

if __name__ == "__main__":

    BNP7 = DistevalLibrary('BNP7/disteval/BNP7.json')
    result = BNP7(parameters={"s": 9., "t": -2.5, "a0": 30./40., "a1": 42./107., "a2": 51./65., "a3": 67./89., "a4": 79./55., "a5": 88./33., "a6": 91./33.})

    print('Numerical Result')
    print(result)
