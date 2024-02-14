import os
import sys
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt

from pyomo.environ import ConcreteModel, Var, Constraint
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.solvers.get_solver import get_solver


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
__author__ = "Kurban Sitterley"

"""
Demonstration file used to make chlorination surrogate model.
"""

surrogate_file = "chlorination_surrogate-test.json"

chlorination_data = (
    os.path.abspath(os.path.join(__location__, os.pardir))
    + "/data/chlorination_cost.csv"
)

solver = get_solver()
data = pd.read_csv(chlorination_data)


def main():
    surrogate = build_surrogate()
    test_surrogate(surrogate)


def build_surrogate(delete_surrogate=True):

    input_labels = ["dose", "flow_mgd"]
    input_bounds = dict(dose=[0.9, 26], flow_mgd=[0.9, 26])
    output_labels = ["capital_cost"]
    data_training, _ = split_training_validation(data, 0.95, seed=len(data))
    stream = StringIO()
    oldstdout = sys.stdout
    sys.stdout = stream
    trainer = PysmoRBFTrainer(
        input_labels=input_labels,
        output_labels=output_labels,
        training_dataframe=data_training,
    )
    trainer.config.basis_function = "gaussian"
    trainer.config.solution_method = "algebraic"
    trainer.config.regularization = True

    trained_rbf = trainer.train_surrogate()
    os.remove("solution.pickle")

    surrogate = PysmoSurrogate(trained_rbf, input_labels, output_labels, input_bounds)
    surrogate.save_to_file(surrogate_file, overwrite=True)
    if delete_surrogate:
        os.remove(surrogate_file)
    sys.stdout = oldstdout
    return surrogate


def test_surrogate(surrogate):
    m = ConcreteModel()
    m.dose = Var(initialize=1)
    m.flow_mgd = Var(initialize=1)
    m.capital_cost = Var(initialize=1e5)

    m.surrogate_blk = SurrogateBlock(concrete=True)
    m.surrogate_blk.build_model(
        surrogate, input_vars=[m.dose, m.flow_mgd], output_vars=[m.capital_cost]
    )

    m.flow_mgd.fix(10)
    m.dose.fix(25)

    _ = solver.solve(m)
    pred_capital = list()

    for dose, flow in zip(data.dose, data.flow_mgd):
        m.flow_mgd.fix(flow)
        m.dose.fix(dose)
        results = solver.solve(m)
        print(f"Termination = {results.solver.termination_condition}")
        pred_capital.append(m.capital_cost.value)

    _, ax = plt.subplots()
    ax.scatter(data.capital_cost, pred_capital, marker=".")
    ax.set_xlabel("Real")
    ax.set_ylabel("Predicted")
    ax.set_title("Chlorination surrogate model - WT3")
    plt.show()


if __name__ == "__main__":
    main()
    plt.close()
