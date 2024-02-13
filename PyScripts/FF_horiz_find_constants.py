import numpy as np
from scipy.optimize import fsolve

from sqlalchemy import create_engine, Column, Integer, Float, MetaData
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker


def bs_integral(Y0, D, sag, d_P, r, zspace):

    """
    This function numerically calculates the biot savart integral for a catenary
    the return is a list of results, one result for each wire
    """

    # define Hfunc which has to be solved to find a and b
    def Hfunc(a):
        return D - 2*a*np.arccosh((sag + a)/a)

    # solve Hfunc to find a and b
    a = fsolve(Hfunc, np.array([100]))[0]
    b = a - (Y0 - sag)

    # now define the function  y = catenary(z)
    def catenary(z):
        return a*np.cosh(z/a) - b

    """
    B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
    """
    ret = 0
    # first calculate the integral for each phase, B will then scale with I as I changes over time
    for p_n, (P_x, I) in enumerate(zip([-d_P, 0, d_P], [np.sin(np.deg2rad(-30)), 0, np.sin(np.deg2rad(-30))])):

        integral_result = np.zeros(3)
        # starting point of curve
        last_P = np.array([P_x, catenary(0), 0])
        for z in zspace[1:]:
            next_P = np.array([P_x, catenary(z), z])
            dl = next_P - last_P
            OP = last_P - r
            OP_length = np.sqrt(np.sum(np.power(OP, 2)))
            OP_unit = OP/OP_length
            integral_result += np.cross(dl, OP_unit)/(OP_length**2)
            last_P = next_P
        # since integral only uses half the wire the other half is same except mirrored over XY plane
        # thus the Z component drops out but X and Y component doubles
        # integral_result = 2*integral_result

        ret += 2*I*integral_result[0]
        # now that the integral has been computed the B can be simulated over time

    return ret


def calculate(Y0, D, sag, d_P):

    zspace = np.linspace(0, D, 200)

    a = bs_integral(Y0, D, sag, d_P, r=np.array([0, 1, 0]), zspace=zspace)
    b = bs_integral(Y0, D, sag, d_P, r=np.array([0, 0.5, 0]), zspace=zspace)
    return b/a

Base = declarative_base()

class CalculationResult(Base):
    __tablename__ = 'calculation_results'
    id = Column(Integer, primary_key=True)
    Y0 = Column(Float)
    D = Column(Float)
    sag = Column(Float)
    d_P = Column(Float)
    result = Column(Float)

class Calculator:
    def __init__(self):
        # SQLite database in memory for demonstration purposes
        self.engine = create_engine('sqlite:///:memory:')
        Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)

    def calculate(self, Y0, D, sag, d_P):
        # Check if the result is already in the database
        session = self.Session()
        result_entry = session.query(CalculationResult).filter_by(
            Y0=Y0, D=D, sag=sag, d_P=d_P
        ).first()

        if result_entry:
            print("Result already computed:", result_entry.result)
        else:
            # Perform the calculation
            result = calculate(Y0=Y0, D=D, sag=sag, d_P=d_P)

            # Store the result in the database
            new_result_entry = CalculationResult(
                Y0=Y0, D=D, sag=sag, d_P=d_P, result=result
            )
            session.add(new_result_entry)
            session.commit()

            print("Calculated result:", result)



# Example usage:
calculator = Calculator()

# First calculation
calculator.calculate(15, 250, 5, 7)

# Different calculation
calculator.calculate(16, 250, 5, 7)

# Repeating the same calculation
calculator.calculate(15, 250, 5, 7)