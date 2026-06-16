import re
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
from scipy.constants import R  # 8.314472  J mol^-1 K^-1

from ...utils import logger

try:
    import numexpr

    HAVE_NUMEXPR = True
except ImportError:
    HAVE_NUMEXPR = False

#from .util import logsumexp

KCAL_MOL_2_KJ_MOL = 4.184
TEMPERATURE = 310
KB = R / 1000
BETA = 1 / (KB * TEMPERATURE)  # mol / KJ

def logsumexp(a, axis=None, b=None, use_numexpr=True):
    """Compute the log of the sum of exponentials of input elements.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : None or int, optional, default=None
        Axis or axes over which the sum is taken. By default `axis` is None,
        and all elements are summed.
    b : array-like, optional
        Scaling factor for exp(`a`) must be of the same shape as `a` or
        broadcastable to `a`.
    use_numexpr : bool, optional, default=True
        If True, use the numexpr library to speed up the calculation, which
        can give a 2-4X speedup when working with large arrays.

    Returns
    -------
    res : ndarray
        The result, ``log(sum(exp(a)))`` calculated in a numerically
        more stable way. If `b` is given then ``log(sum(b*exp(a)))``
        is returned.

    See Also
    --------
    numpy.logaddexp, numpy.logaddexp2, scipy.misc.logsumexp (soon to be replaced with  scipy.special.logsumexp)

    Notes
    -----
    This is based on scipy.misc.logsumexp but with optional numexpr
    support for improved performance.
    """

    a = np.asarray(a)

    a_max = np.amax(a, axis=axis, keepdims=True)

    if a_max.ndim > 0:
        a_max[~np.isfinite(a_max)] = 0
    elif not np.isfinite(a_max):
        a_max = 0

    if b is not None:
        b = np.asarray(b)
        if use_numexpr and HAVE_NUMEXPR:
            out = np.log(numexpr.evaluate("b * exp(a - a_max)").sum(axis))
        else:
            out = np.log(np.sum(b * np.exp(a - a_max), axis=axis))
    else:
        if use_numexpr and HAVE_NUMEXPR:
            out = np.log(numexpr.evaluate("exp(a - a_max)").sum(axis))
        else:
            out = np.log(np.sum(np.exp(a - a_max), axis=axis))

    a_max = np.squeeze(a_max, axis=axis)
    out += a_max

    return out



class CSVParser:
    def __init__(self, csv, temperature: float = 310):
        self.csv = Path(csv)
        self.BETA = 1 / (KB * temperature)

    def parse_df(self) -> pd.DataFrame:
        df = pd.read_csv(self.csv, na_filter=True, memory_map=True, sep=r"\s+")

        # remove Time, U, pV column, only save the dU column
        u_k = df.loc[:, ~df.columns.isin(["Time", "U", "pV"])]

        # add pV and convert to reduced potential
        u_k = u_k.add(df.pV, axis="rows") * self.BETA

        # rename window list
        u_k.columns = [f"{i}" for i in range(len(u_k.columns))]

        # add a column to df for using groupby
        u_k["lambda"] = f"lambda_{self.csv.parent.name}"
        u_k["window"] = f"{int(self.csv.parent.name)}"  # for csv
        # set lambda index for later groupby
        return u_k.set_index(["lambda", "window"])

class XVGParser:
    def __init__(self, xvg: Union[str, Path]) -> None:
        self.xvg = Path(xvg)

    def _extract_header(self):
        self.headlines = []
        self.header_cnt = 0
        with open(self.xvg) as f:
            for line in f:
                if line.startswith("#"):
                    self.header_cnt += 1
                    continue
                elif line.startswith("@"):
                    self.headlines.append(line)
                else:
                    break
                self.header_cnt += 1

    def parse_header(self):
        self._extract_header()
        re_T = re.compile(r"T\s+=\s+(\d+)")
        re_state = re.compile(r"state\s+(\d+)")
        re_lambda_1 = re.compile(r"\(([^\(]*?)\)\s+=\s+\(([^\(]*?)\)")
        re_lambda_2 = re.compile(r"(vdw-lambda\s+)=\s+([0-9.]+)")
        re_s_legend = re.compile(r's\d+ legend\s+"(.*)"')
        self.columns = []
        for line in self.headlines:
            if "subtitle" in line:
                # extract temperature
                t_matched = re_T.search(line)
                state_matched = re_state.search(line)
                lambda_matched = re_lambda_1.search(line)
                lambda_2_matched = re_lambda_2.search(line)
                if t_matched:
                    self.temperature = float(t_matched.group(1))
                else:
                    self.temperature = TEMPERATURE
                if state_matched:
                    self.state = int(state_matched.group(1))
                else:
                    raise RuntimeError("state cannot find in the xvg file, the molecule_dynamics has failed")
                if lambda_matched:
                    self.lambda_key = lambda_matched.group(1).split(",")
                    self.lambda_value = list(map(float, lambda_matched.group(2).split(",")))
                elif lambda_2_matched:
                    self.lambda_key = [lambda_2_matched.group(1)]
                    self.lambda_value = [float(lambda_2_matched.group(2))]
                else:
                    raise RuntimeError("Cannot not find the lambda info in xvg file")
            s_match = re_s_legend.search(line)
            if s_match:
                self.columns.append(s_match.group(1))

    def _extract_data(self):
        self.parse_header()
        cols = ["time"] + [
            col + "{}[duplicated]".format(i) if col in self.columns[:i] else col for i, col, in enumerate(self.columns)
        ]
        df = pd.read_csv(
            self.xvg,
            sep=r"\s+",
            header=None,
            skiprows=self.header_cnt,
            na_filter=True,
            memory_map=True,
            names=cols,
            float_precision="high",
        )
        df = df[df.columns[~df.columns.str.endswith("[duplicated]")]]
        return df

    def extract_dh(self) -> pd.DataFrame:
        df = self._extract_data()
        h_col_match = r"\xD\f{}H \xl\f{}"
        pv_col_match = "pV"

        # grab only dH columns
        DHcols = [col for col in df.columns if (h_col_match in col)]
        dH = df[DHcols]

        # gromacs also gives us pV directly; need this for reduced potential
        pv_cols = [col for col in df.columns if (pv_col_match in col)]

        u_k = dict()
        cols = []

        for col in dH:
            u_col = eval(col.split("to")[1])
            # calculate reduced potential u_k = dH + pV + U
            u_k[u_col] = BETA * dH[col].values
            if pv_cols:
                u_k[u_col] += BETA * df[pv_cols[0]].values
            cols.append(u_col)
        u_k = pd.DataFrame(u_k, columns=cols)
        u_k["time"] = df["time"]
        for i, l in enumerate(self.lambda_key):
            u_k[l] = self.lambda_value[i]
        newidx = ["time"] + self.lambda_key
        return u_k.set_index(newidx)

    def extract_dhdl(self) -> pd.DataFrame:
        dh_dl_col_match = r"dH/d\xl\f{}"
        df = self._extract_data()
        # grab only lambda columns
        lambda_cols = [col for col in df.columns if dh_dl_col_match in col]
        dHdl = df[lambda_cols]

        u_k = dict()
        cols = []

        for col, lambda_name in zip(dHdl, self.lambda_key):
            newcol = lambda_name.split("-")[0]
            u_k[newcol] = BETA * dHdl[col]
            cols.append(newcol)
        u_k = pd.DataFrame(u_k, columns=cols)
        u_k["time"] = df["time"]
        for i, l in enumerate(self.lambda_key):
            u_k[l] = self.lambda_value[i]
        newidx = ["time"] + self.lambda_key
        return u_k.set_index(newidx)


class BarForDDG:
    def __init__(self, u_nks: List[pd.DataFrame], number_blocks=5) -> None:
        """
        u_nks contains a list of dataframes which has a length of n (n = number_of_lambda_window)
        in each dataframe, the reducued potential of du must be specfied
        such as:
                      dU_0  dU_1       dU_2       dU_3       dU_4 ...
        lambda
        lambda_1 -5.468750   0.0  14.886719  28.898438  42.105469 ...
        lambda_1 -7.328125   0.0  17.675781  34.675781  51.242188 ...
        lambda_1 -6.562500   0.0  15.859375  30.917969  45.667969 ...
        lambda_1 -8.828125   0.0  20.238281  39.328125  57.757812 ...
        lambda_1 -5.808594   0.0  15.292969  29.664062  43.164062 ...
        ...            ...   ...        ...        ...        ... ...
        lambda_1 -5.671875   0.0  15.410156  30.003906  43.863281 ...
        lambda_1 -8.082031   0.0  19.316406  37.605469  55.199219 ...
        lambda_1 -7.828125   0.0  19.066406  37.449219  55.359375 ...
        lambda_1 -6.242188   0.0  16.320312  32.050781  47.386719 ...
        lambda_1 -8.089844   0.0  18.871094  36.750000  53.578125 ...

        the lambda column is used for the index and it is used for later groupby function.
        """
        self.number_blocks = number_blocks
        # df has all of the data
        df = pd.concat(u_nks)
        self.number_lambda = df.shape[1]
        self.number_frame = df.shape[0] / df.shape[1]
        self.lambda_name = df.columns.values.tolist()
        self.groups = df.groupby(level=df.index.names[1])

    @classmethod
    def groupby_xvg(cls, xvg_path, xvg_name="prod_npt.xvg"):
        xvgs = sorted(Path(xvg_path).glob(f"*/{xvg_name}"), key=lambda x: int(x.parent.name))
        u_nks = [XVGParser(xvg).extract_dh() for xvg in xvgs]

    @classmethod
    def groupby_csv(cls, csv_path, csv_name="prod_npt.csv"):
        csvs = sorted(Path(csv_path).glob(f"*/{csv_name}"), key=lambda x: int(x.parent.name))
        u_nks = [CSVParser(csv).parse_df() for csv in csvs]
        return cls(u_nks)

    def run(self):
        dF, dF_uncertainty = [], []
        sa, sb = [], []
        for k in range(self.number_lambda - 1):
            sub_df = self.groups.get_group(self.lambda_name[k])
            wf = sub_df.iloc[:, k + 1] - sub_df.iloc[:, k]
            sub_df2 = self.groups.get_group(self.lambda_name[k + 1])
            wr = sub_df2.iloc[:, k] - sub_df2.iloc[:, k + 1]
            result = self._solve_bar_equation(wf, wr)
            dF.append(result[0])
            dF_uncertainty.append(result[1])
            sa.append((np.sum(wf) / len(wf) - result[0]))
            sb.append((np.sum(wr) / len(wr) + result[0]))

        # five blocks to compute the std deviation
        dFs = []
        stride = self.number_frame // self.number_blocks
        for i in range(self.number_blocks):
            start = int(stride * i)
            end = int(stride * (i + 1))
            dF_block = []
            for k in range(self.number_lambda - 1):
                sub_df = self.groups.get_group(self.lambda_name[k])
                # wf = sub_df.iloc[i::n_blocks, k + 1] - sub_df.iloc[i::n_blocks, k]
                wf = sub_df.iloc[start:end, k + 1] - sub_df.iloc[start:end, k]
                sub_df2 = self.groups.get_group(self.lambda_name[k + 1])
                # wr = sub_df2.iloc[i::n_blocks, k] - sub_df2.iloc[i::n_blocks, k + 1]
                wr = sub_df2.iloc[start:end, k] - sub_df2.iloc[start:end, k + 1]
                result = self._solve_bar_equation(wf, wr, compute_uncertainty=True)
                dF_block.append(result[0])
            dFs.append(dF_block)
        return (
            np.array(dF) / BETA / KCAL_MOL_2_KJ_MOL,
            np.array(dF_uncertainty) / BETA / KCAL_MOL_2_KJ_MOL,
            np.std(dFs, axis=0, ddof=1) / BETA / KCAL_MOL_2_KJ_MOL,
            sa,
            sb,
        )
        # np.sum(dF) / BETA / KCAL_MOL_2_KJ_MOL,
        # max(np.sqrt(np.sum(np.std(dFs, axis=0) ** 2)) / BETA / KCAL_MOL_2_KJ_MOL,
        #     df_sum_uncertainty)

    def _solve_bar_equation(
        self,
        w_F,
        w_R,
        DeltaF=0.0,
        compute_uncertainty=True,
        maximum_iterations=500,
        relative_tolerance=1.0e-12,
        verbose=False,
        return_dict=False,
    ) -> Tuple[float, float]:
        """Compute free energy difference using the Bennett acceptance ratio (BAR) method.

        Parameters
        ----------
        w_F : np.ndarray
            w_F[t] is the forward work value from snapshot t.
            t = 0...(T_F-1)  Length T_F is deduced from vector.
        w_R : np.ndarray
            w_R[t] is the reverse work value from snapshot t.
            t = 0...(T_R-1)  Length T_R is deduced from vector.
        DeltaF : float, optional, default=0.0
            DeltaF can be set to initialize the free energy difference with a guess
        compute_uncertainty : bool, optional, default=True
            if False, only the free energy is returned
        uncertainty_method: string, optional, default=BAR
            There are two possible uncertainty estimates for BAR.  One agrees with MBAR for two states exactly;
            The other only agrees with MBAR in the limit of good overlap. See below.
        maximum_iterations : int, optional, default=500
            can be set to limit the maximum number of iterations performed
        relative_tolerance : float, optional, default=1E-11
            can be set to determine the relative tolerance convergence criteria (defailt 1.0e-11)
        verbose : bool
            should be set to True if verbse debug output is desired (default False)
        method : str, optional, defualt='false-position'
            choice of method to solve BAR nonlinear equations, one of 'self-consistent-iteration' or 'false-position' (default: 'false-position')
        iterated_solution : bool, optional, default=True
            whether to fully solve the optimized BAR equation to consistency, or to stop after one step, to be
            equivalent to transition matrix sampling.
        return_dict : bool, default False
            If true, returns are a dict, else they are a tuple

        Returns
        -------
        'Delta_f' : float
            Free energy difference
            If return_dict, key is 'Delta_f'
        'dDelta_f': float
            Estimated standard deviation of free energy difference
            If return_dict, key is 'dDelta_f'


        References
        ----------

        [1] Shirts MR, Bair E, Hooker G, and Pande VS. Equilibrium free energies from nonequilibrium
        measurements using maximum-likelihood methods. PRL 91(14):140601, 2003.

        Notes
        -----
        The false position method is used to solve the implicit equation.

        Examples
        --------
        Compute free energy difference between two specified samples of work values.

        >>> from pymbar import testsystems
        >>> [w_F, w_R] = testsystems.gaussian_work_example(mu_F=None, DeltaF=1.0, seed=0)
        >>> results = BAR(w_F, w_R, return_dict=True)
        >>> print('Free energy difference is {:.3f} +- {:.3f} kT'.format(results['Delta_f'], results['dDelta_f']))
        Free energy difference is 1.088 +- 0.050 kT

        Test completion of various other schemes.

        >>> results = BAR(w_F, w_R, method='self-consistent-iteration', return_dict=True)
        >>> results = BAR(w_F, w_R, method='false-position', return_dict=True)
        >>> results = BAR(w_F, w_R, method='bisection', return_dict=True)

        """

        result_vals = dict()

        UpperB = -(logsumexp(-w_F) - np.log(w_F.size))
        LowerB = logsumexp(-w_R) - np.log(w_R.size)

        FUpperB = self._BARzero(w_F, w_R, UpperB)
        FLowerB = self._BARzero(w_F, w_R, LowerB)
        nfunc = 2

        if np.isnan(FUpperB) or np.isnan(FLowerB):
            # this data set is returning NAN -- will likely not work.  Return 0, print a warning:
            # consider returning more information about failure
            logger.warning(
                "BAR is likely to be inaccurate because of poor overlap. "
                "Improve the sampling, or decrease the spacing betweeen states."
                "For now, guessing that the free energy difference is 0 with no uncertainty."
            )
            if compute_uncertainty:
                result_vals["Delta_f"] = 0.0
                result_vals["dDelta_f"] = 0.0
                if return_dict:
                    return result_vals
                return 0.0, 0.0
            else:
                result_vals["Delta_f"] = 0.0
                if return_dict:
                    return result_vals
                return 0.0

        while FUpperB * FLowerB > 0:
            # if they have the same sign, they do not bracket.  Widen the bracket until they have opposite signs.
            # There may be a better way to do this, and the above bracket should rarely fail.
            FAve = (UpperB + LowerB) / 2
            UpperB = UpperB - max(abs(UpperB - FAve), 0.1)
            LowerB = LowerB + max(abs(LowerB - FAve), 0.1)
            FUpperB = self._BARzero(w_F, w_R, UpperB)
            FLowerB = self._BARzero(w_F, w_R, LowerB)
            nfunc += 2

        # Iterate to convergence or until maximum number of iterations has been exceeded.

        for iteration in range(maximum_iterations):

            DeltaF_old = DeltaF
            # Predict the new value
            if (LowerB == 0.0) and (UpperB == 0.0):
                DeltaF = 0.0
                FNew = 0.0
            else:
                DeltaF = UpperB - FUpperB * (UpperB - LowerB) / (FUpperB - FLowerB)
                FNew = self._BARzero(w_F, w_R, DeltaF)
                nfunc += 1

            if FNew == 0:
                # Convergence is achieved.
                if verbose:
                    logger.info("Convergence achieved.")
                relative_change = 10 ** (-15)
                break

            # Check for convergence.
            if DeltaF == 0.0:
                # The free energy difference appears to be zero -- return.
                if verbose:
                    logger.info("The free energy difference appears to be zero.")
                break

            relative_change = abs((DeltaF - DeltaF_old) / DeltaF)
            if verbose:
                logger.info("relative_change = {:12.3f}".format(relative_change))

            if (iteration > 0) and (relative_change < relative_tolerance):
                # Convergence is achieved.
                if verbose:
                    logger.info("Convergence achieved.")
                break
            if FUpperB * FNew < 0:
                # these two now bracket the root
                LowerB = DeltaF
                FLowerB = FNew
            elif FLowerB * FNew <= 0:
                # these two now bracket the root
                UpperB = DeltaF
                FUpperB = FNew
            else:
                logger.warning("WARNING: Cannot determine bound on free energy")

            if verbose:
                logger.info("iteration {:5d}: DeltaF = {:16.3f}".format(iteration, DeltaF))

        # Report convergence, or warn user if not achieved.
        if iteration < maximum_iterations:
            if verbose:
                logger.info(
                    "Converged to tolerance of {:e} in {:d} iterations ({:d} function evaluations)".format(
                        relative_change, iteration, nfunc
                    )
                )
        else:
            logger.warning(
                "Did not converge to within specified tolerance. max_delta = {:f}, TOLERANCE = {:f}, MAX_ITS = %d".format(
                    relative_change, relative_tolerance, maximum_iterations
                )
            )

        """
        Compute asymptotic variance estimate using Eq. 10a of Bennett,
        1976 (except with n_1<f>_1^2 in the second denominator, it is
        an error in the original NOTE: The 'BAR' and 'MBAR' estimators
        do not agree for poor overlap. This is not because of
        numerical precision, but because they are fundamentally
        different estimators. For poor overlap, 'MBAR' diverges high,
        and 'BAR' diverges by being too low. In situations they are
        noticeably from each other, they are also pretty different
        from the true answer (obtained by calculating the standard
        deviation over lots of realizations).

        First, we examine the 'BAR' equation. Rederive from Bennett, substituting (8) into (7)

        (8)    -> W = [q0/n0 exp(-U1) + q1/n1 exp(-U0)]^-1
                    <(W exp(-U1))^2 >_0         <(W exp(-U0))^2 >_1
        (7)    -> -----------------------  +   -----------------------   - 1/n0 - 1/n1
                   n_0 [<(W exp(-U1)>_0]^2      n_1 [<(W exp(-U0)>_1]^2

            Const cancels out of top and bottom.   Wexp(-U0) = [q0/n0 exp(-(U1-U0)) + q1/n1]^-1
                                                             =  n1/q1 [n1/n0 q0/q1 exp(-(U1-U0)) + 1]^-1
                                                             =  n1/q1 [exp (M+(F1-F0)-(U1-U0)+1)^-1]
                                                             =  n1/q1 f(x)
                                                   Wexp(-U1) = [q0/n0 + q1/n1 exp(-(U0-U1))]^-1
                                                             =  n0/q0 [1 + n0/n1 q1/q0 exp(-(U0-U1))]^-1
                                                             =  n0/q0 [1 + exp(-M+[F0-F1)-(U0-U1))]^-1
                                                             =  n0/q0 f(-x)


                  <(W exp(-U1))^2 >_0          <(W exp(-U0))^2 >_1
         (7) -> -----------------------   +  -----------------------   - 1/n0 - 1/n1
                n_0 [<(W exp(-U1)>_0]^2      n_1 [<(W exp(-U0)>_1]^2

                   <[n0/q0 f(-x)]^2>_0        <[n1/q1 f(x)]^2>_1
                -----------------------  +  ------------------------   -1/n0 -1/n1
                  n_0 <n0/q0 f(-x)>_0^2      n_1 <n1/q1 f(x)>_1^2

               1      <[f(-x)]^2>_0                 1        <[f(x)]^2>_1
               -  [-----------------------  - 1]  + -  [------------------------  - 1]
               n0      <f(-x)>_0^2                  n1      n_1<f(x)>_1^2

        where f = the fermi function, 1/(1+exp(-x))

        This formula the 'BAR' equation works for works for free
        energies (F0-F1) that don't satisfy the BAR equation.  The
        'MBAR' equation, detailed below, only works for free energies
        that satisfy the equation.


        Now, let's look at the MBAR version of the uncertainty.  This
        is written (from Shirts and Chodera, JPC, 129, 124105, Equation E9) as

              [ n0<f(x)f(-x)>_0 + n1<f(x)f(-x)_1 ]^-1 - n0^-1 - n1^-1

              we note the f(-x) + f(x)  = 1, and change this to:

              [ n0<(1-f(-x)f(-x)>_0 + n1<f(x)(1-f(x))_1 ]^-1 - n0^-1 - n1^-1

              [ n0<f(-x)-f(-x)^2)>_0 + n1<f(x)-f(x)^2)_1 ]^-1 - n0^-1 - n1^-1

                                                1                                         1     1
              --------------------------------------------------------------------    -  --- - ---
                 n0 <f(-x)>_0 - n0 <[f(-x)]^2>_0 + n1 <f(x)>_1 + n1 <[f(x)]^2>_1          n0    n1


        Removing the factor of - (T_F + T_R)/(T_F*T_R)) from both, we compare:

                  <[f(-x)]^2>_0          <[f(x)]^2>_1
              [------------------]  + [---------------]
                 n0 <f(-x)>_0^2          n1 <f(x)>_1^2

                                                1
              --------------------------------------------------------------------
                 n0 <f(-x)>_0 - n0 <[f(-x)]^2>_0 + n1 <f(x)>_1 + n1 <[f(x)]^2>_1

        denote: <f(-x)>_0 = afF
                <f(-x)^2>_0 = afF2
                <f(x)>_1 = afR
                <f(x)^2>_1 = afF2

        Then we can look at both of these as:

        variance_BAR = (afF2/afF**2)/T_F + (afR2/afR**2)/T_R
        variance_MBAR = 1/(afF*T_F - afF2*T_F + afR*T_R - afR2*T_R)

        Rearranging:

        variance_BAR = (afF2/afF**2)/T_F + (afR2/afR**2)/T_R
        variance_MBAR = 1/(afF*T_F + afR*T_R - (afF2*T_F +  afR2*T_R))

        # check the steps below?  Not quite sure.
        variance_BAR = (afF2/afF**2) + (afR2/afR**2)  = (afF2 + afR2)/afR**2
        variance_MBAR = 1/(afF + afR - (afF2 +  afR2)) = 1/(2*afR-(afF2+afR2))

        Definitely not the same.  Now, the reason that they both work
        for high overlap is still not clear. We will determine the
        difference at some point.

        see https://github.com/choderalabpymbar/issues/281 for more information.

        Now implement the two computations.
        """

        if compute_uncertainty:
            T_F = float(w_F.size)  # number of forward work values
            T_R = float(w_R.size)  # number of reverse work values
            M = np.log(T_F / T_R)
            C = M - DeltaF

            sum_cosh = np.sum(1.0 / (2 + 2 * np.cosh(w_F + C))) + np.sum(1.0 / (2 + 2 * np.cosh(-w_R + C)))
            dDeltaF = np.sqrt((1.0 / (sum_cosh / (T_F + T_R)) - (T_F + T_R) / T_F - (T_F + T_R) / T_R) / (T_F + T_R))
            return DeltaF, dDeltaF
        else:
            return DeltaF, 0

        # if compute_uncertainty:
        #     # Determine number of forward and reverse work values provided.
        #     T_F = float(w_F.size)  # number of forward work values
        #     T_R = float(w_R.size)  # number of reverse work values

        #     # Compute log ratio of forward and reverse counts.
        #     M = np.log(T_F / T_R)

        #     C = M - DeltaF

        #     # In theory, overflow handling should not be needed now, because we use numlogexp or a custom routine?

        #     # fF = 1 / (1 + np.exp(w_F + C)), but we need to handle overflows
        #     exp_arg_F = (w_F + C)
        #     max_arg_F  = np.max(exp_arg_F)
        #     log_fF = - np.log(np.exp(-max_arg_F) + np.exp(exp_arg_F - max_arg_F))
        #     afF  = np.exp(logsumexp(log_fF)-max_arg_F)/T_F

        #     # fR = 1 / (1 + np.exp(w_R - C)), but we need to handle overflows
        #     exp_arg_R = (w_R - C)
        #     max_arg_R  = np.max(exp_arg_R)
        #     log_fR = - np.log(np.exp(-max_arg_R) + np.exp(exp_arg_R - max_arg_R))
        #     afR = np.exp(logsumexp(log_fR)-max_arg_R)/T_R

        #     afF2 = np.exp(logsumexp(2*log_fF)-2*max_arg_F)/T_F
        #     afR2 = np.exp(logsumexp(2*log_fR)-2*max_arg_R)/T_R

        #     nrat = (T_F + T_R)/(T_F * T_R)   # same for both methods

        #     variance = (afF2/afF**2)/T_F + (afR2/afR**2)/T_R - nrat
        #     dDeltaF = np.sqrt(variance)

    def _BARzero(self, w_F, w_R, DeltaF):
        """A function that when zeroed is equivalent to the solution of
        the Bennett acceptance ratio.

        from http://journals.aps.org/prl/pdf/10.1103/PhysRevLett.91.140601
        D_F = M + w_F - Delta F
        D_R = M + w_R - Delta F

        we want:
        \sum_N_F (1+exp(D_F))^-1 = \sum N_R N_R <(1+exp(-D_R))^-1>
        ln \sum N_F (1+exp(D_F))^-1>_F = \ln \sum N_R exp((1+exp(-D_R))^(-1)>_R
        ln \sum N_F (1+exp(D_F))^-1>_F - \ln \sum N_R exp((1+exp(-D_R))^(-1)>_R = 0

        Parameters
        ----------
        w_F : np.ndarray
            w_F[t] is the forward work value from snapshot t.
            t = 0...(T_F-1)  Length T_F is deduced from vector.
        w_R : np.ndarray
            w_R[t] is the reverse work value from snapshot t.
            t = 0...(T_R-1)  Length T_R is deduced from vector.
        DeltaF : float
            Our current guess

        Returns
        -------
        fzero : float
            a variable that is zeroed when DeltaF satisfies BAR.

        Examples
        --------
        Compute free energy difference between two specified samples of work values.

        >>> from pymbar import testsystems
        >>> [w_F, w_R] = testsystems.gaussian_work_example(mu_F=None, DeltaF=1.0, seed=0)
        >>> DeltaF = _BARzero(w_F, w_R, 0.0)

        """

        np.seterr(over="raise")  # raise exceptions to overflows
        w_F = np.array(w_F, np.float64)
        w_R = np.array(w_R, np.float64)
        DeltaF = float(DeltaF)

        # Recommended stable implementation of BAR.

        # Determine number of forward and reverse work values provided.
        T_F = float(w_F.size)  # number of forward work values
        T_R = float(w_R.size)  # number of reverse work values

        # Compute log ratio of forward and reverse counts.
        M = np.log(T_F / T_R)

        # Compute log numerator. We have to watch out for overflows.  We
        # do this by making sure that 1+exp(x) doesn't overflow, choosing
        # to always exponentiate a negative number.

        # log f(W) = - log [1 + exp((M + W - DeltaF))]
        #          = - log ( exp[+maxarg] [exp[-maxarg] + exp[(M + W - DeltaF) - maxarg]] )
        #          = - maxarg - log(exp[-maxarg] + exp[(M + W - DeltaF) - maxarg])
        # where maxarg = max((M + W - DeltaF), 0)

        exp_arg_F = M + w_F - DeltaF
        # use boolean logic to zero out the ones that are less than 0, but not if greater than zero.
        max_arg_F = np.choose(np.less(0.0, exp_arg_F), (0.0, exp_arg_F))
        try:
            log_f_F = -max_arg_F - np.log(np.exp(-max_arg_F) + np.exp(exp_arg_F - max_arg_F))
        except:
            # give up; if there's overflow, return zero
            logger.error("The input data results in overflow in BAR")
            return np.nan
        log_numer = logsumexp(log_f_F)

        # Compute log_denominator.
        # log f(R) = - log [1 + exp(-(M + W - DeltaF))]
        #          = - log ( exp[+maxarg] [exp[-maxarg] + exp[(M + W - DeltaF) - maxarg]] )
        #          = - maxarg - log[exp[-maxarg] + (T_F/T_R) exp[(M + W - DeltaF) - maxarg]]
        # where maxarg = max( -(M + W - DeltaF), 0)

        exp_arg_R = -(M - w_R - DeltaF)
        # use boolean logic to zero out the ones that are less than 0, but not if greater than zero.
        max_arg_R = np.choose(np.less(0.0, exp_arg_R), (0.0, exp_arg_R))
        try:
            log_f_R = -max_arg_R - np.log(np.exp(-max_arg_R) + np.exp(exp_arg_R - max_arg_R))
        except:
            logger.error("The input data results in overflow in BAR")
            return np.nan
        log_denom = logsumexp(log_f_R)

        # This function must be zeroed to find a root
        fzero = log_numer - log_denom

        # return options to standard settings so we don't disturb other functionality.
        np.seterr(over="warn")
        return fzero


if __name__ == "__main__":
    from pathlib import Path

    # p = Path("/home/haomiao/work/fep/fep_analyze/test/MCL1/27_to_45/rbfe/stage1/")
    # bar = BarEstimator.groupby_xvg(p)
    # d1, d2, d3 = bar.analyze()
    # print("xvg file bar")
    # print(d1)
    # bar = BarEstimator.groupby_csv(p)
    # d1, d2, d3 = bar.analyze()
    # print("csv file bar")
    # print(d1)

    p = Path("/home/haomiao/temp/bar/stage1")
    bar = BarForDDG.groupby_csv(p)
    d, _, _ = bar.run()
    print(d)
