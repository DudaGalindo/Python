from scipy.linalg import eig
import numpy as np

def eigV(A, B=None):
        if B is None:
            evalues, evectors = eig(A)
        else:
            evalues, evectors = eig(A, B)

        if all(eigs == 0 for eigs in evalues.imag):
            if all(eigs > 0 for eigs in evalues.real):
                idxp = evalues.real.argsort()  # positive in increasing order
                idxn = np.array([], dtype=int)
            else:
                # positive in increasing order
                idxp = evalues.real.argsort()[int(len(evalues) / 2):]
                # negative in decreasing order
                idxn = evalues.real.argsort()[int(len(evalues) / 2) - 1:: -1]

        else:
            # positive in increasing order
            idxp = evalues.imag.argsort()[int(len(evalues) / 2):]
            # negative in decreasing order
            idxn = evalues.imag.argsort()[int(len(evalues) / 2) - 1:: -1]

        idx = np.hstack([idxp, idxn])

        return evalues[idx], evectors[:, idx]
