from __future__ import annotations
from dataclasses import dataclass
from typing import Any, List, Literal, Union
import operator
import pandas as pd
import numpy as np
import re

__all__ = ['select_phase']

# define allowed operations
OPERATIONS = Literal['==', '!=', '<', '<=', '>', '>=', 'in', 'not in']
OPERATION_MAPPER = {
    '==': operator.eq,
    '!=': operator.ne,
    '<': operator.lt,
    '<=': operator.le,
    '>': operator.gt,
    '>=': operator.ge,
}

# define data structures for conditions and groups of conditions
@dataclass
class Condition:
    col: str
    op: OPERATIONS
    value: Any = None
    negate: bool = False

@dataclass
class ConditionGroup:
    logic: Literal['and', 'or']
    conditions: List[Union['Condition', 'ConditionGroup']]

@dataclass
class Not:
    group: ConditionGroup


# predefined phases based on common astrophysical events
DEFAULT_PHASES = {
    'BNS_merger': "kstar_1 == 13 & kstar_2 == 13 & evol_type == 3",
    'NSBH_merger': "(kstar_1 == 14 & kstar_2 == 13 & evol_type == 3) | (kstar_1 == 13 & kstar_2 == 14 & evol_type == 3)",
    'BBH_merger': "kstar_1 == 14 & kstar_2 == 14 & evol_type == 3",
    'BH_MS': "(kstar_1 == 14 & kstar_2 in [0,1] & sep > 0) | (kstar_1 in [0,1] & kstar_2 == 14 & sep > 0)",
    'NS_MS': "(kstar_1 == 13 & kstar_2 in [0,1] & sep > 0) | (kstar_1 in [0,1] & kstar_2 == 13 & sep > 0)",
    'WD_MS': "(kstar_1 in [10,11,12] & kstar_2 in [0,1] & sep > 0) | (kstar_1 in [0,1] & kstar_2 in [10,11,12] & sep > 0)",
    'BH_GS': "(kstar_1 == 14 & kstar_2 == 3 & sep > 0) | (kstar_1 == 3 & kstar_2 == 14 & sep > 0)",
    'NS_GS': "(kstar_1 == 13 & kstar_2 == 3 & sep > 0) | (kstar_1 == 3 & kstar_2 == 13 & sep > 0)",
    'WD_GS': "(kstar_1 in [10,11,12] & kstar_2 == 3 & sep > 0) | (kstar_1 == 3 & kstar_2 in [10,11,12] & sep > 0)",
}


def _condition_to_mask(df: pd.DataFrame, condition: Condition) -> pd.Series:
    """Convert a Condition to a boolean mask for the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to apply the condition to
    condition : Condition
        Condition to convert to a mask

    Returns
    -------
    pd.Series
        Boolean mask where the condition is met

    Raises
    ------
    ValueError
        If the operation is not supported
    """
    series = df[condition.col]

    # basic operations are easy, add extra stuff for "in" and "not in"
    if condition.op in OPERATION_MAPPER:
        mask = OPERATION_MAPPER[condition.op](series, condition.value)
    elif condition.op == 'in':
        mask = series.isin(condition.value)
    elif condition.op == 'not in':
        mask = ~series.isin(condition.value)
    else:
        raise ValueError(f"Unsupported operation: {condition.op}")
    
    # flip the condition if negate is set
    if condition.negate:
        mask = ~mask

    return mask

def _group_to_mask(df: pd.DataFrame, group: Union[Condition, ConditionGroup]) -> pd.Series:
    """Convert a ConditionGroup to a boolean mask for the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to apply the group of conditions to
    group : ConditionGroup or Condition
        Group of conditions to convert to a mask

    Returns
    -------
    mask : pd.Series
        Boolean mask where the group of conditions is met

    Raises
    ------
    ValueError
        If the logic is not supported (e.g. not 'and' or 'or')
    """
    if isinstance(group, Condition):
        return _condition_to_mask(df, group)
    
    if group.logic == 'and':
        combined_mask = [_condition_to_mask(df, cond) if isinstance(cond, Condition)
                         else _group_to_mask(df, cond) for cond in group.conditions]
        return np.logical_and.reduce(combined_mask)
    elif group.logic == 'or':
        combined_mask = [_condition_to_mask(df, cond) if isinstance(cond, Condition)
                         else _group_to_mask(df, cond) for cond in group.conditions]
        return np.logical_or.reduce(combined_mask)
    else:
        raise ValueError(f"Unsupported logic: {group.logic}")
    
def _parse_condition(cond: Union[str, Condition]) -> Condition:
    """Parse a single condition string into a Condition object.

    Examples:
        "mass_1 > 4"
        "kstar_2 in [13,14]"
        "sep != 0"

    Parameters
    ----------
    cond : str or Condition
        Condition string to parse or already a Condition object
    Returns
    -------
    parsed_cond : Condition
        Parsed Condition object
    Raises
    ------
    ValueError
        If the condition string is invalid
    """
    if isinstance(cond, Condition):
        return cond

    # all conditions should have at least 3 parts: col, op, value
    parts = cond.split()
    if len(parts) < 3:
        raise ValueError(f"Invalid condition string: {cond}")

    # column names are arbitrary strings, so just take the first part
    col = parts[0].strip()

    # operation must be one of the defined operations, ensure matches OPERATIONS literal
    op = parts[1].strip()
    if op not in OPERATIONS.__args__:
        raise ValueError(f"Invalid operation: {op}")

    # value can be a single number or a list of numbers (for 'in'/'not in')
    value_str = ' '.join(parts[2:])
    if ',' in value_str:
        value_str = value_str.strip('[]() ')
        value = [v.strip() for v in value_str.split(',')]
        assert all(v.replace('.','',1).isdigit() for v in value if v), "All values must be numeric"
        value = [float(v) if '.' in v else int(v) for v in value if v]
    else:
        value = value_str.strip()
        if value.replace('.','',1).isdigit():
            value = float(value) if '.' in value else int(value)
        else:
            raise ValueError("Value must be numeric or a comma-separated list of numerics")

    return Condition(col=col, op=op, value=value)

def __tokenize(expr: str) -> List[str]:
    """Tokenize by parentheses, &, |, ~, and spaces."""
    tokens = re.findall(r'\(|\)|&|\||~|[^&|()~]+', expr)
    return [t.strip() for t in tokens if t.strip()]

def __apply_op(ops, output):
    """Helper function to apply an operator from the ops stack to the output stack."""
    op = ops.pop()
    if op == "~":
        a = output.pop()
        if isinstance(a, ConditionGroup):
            output.append(Not(a))
        else:
            a.negate = not a.negate
            output.append(a)
    else:
        b, a = output.pop(), output.pop()
        logic = "and" if op == "&" else "or"
        output.append(ConditionGroup(logic, [a, b]))
    return ops, output

def _parse_conditions(conds: str) -> Union[ConditionGroup, Condition]:
    """Parse a string of conditions into a ConditionGroup object.
    
    Examples:
        "mass_1 > 4 & mass_2 < 3"
        "(mass_1 > 4 | mass_2 < 3) & evol_type == 3"
        "~(sep == 0) & (kstar_1 == 14 | kstar_2 == 14)"

    Parameters
    ----------
    conds : str
        String of conditions to parse

    Returns
    -------
    parsed_group : ConditionGroup or Condition
        Parsed ConditionGroup object

    Raises
    ------
    ValueError
        If the condition string is invalid
    """
    tokens = __tokenize(conds)
    precedence = {"~":3, "&":2, "|":1}
    output, ops = [], []

    # loop over each token
    for token in tokens:
        if token == "(":
            ops.append(token)
        elif token == ")":
            # when you reach the end of a group, apply all ops until "("
            while ops and ops[-1] != "(":
                ops, output = __apply_op(ops, output)
            if not ops:
                raise ValueError("Mismatched parentheses")
            ops.pop()
        elif token in ("&", "|", "~"):
            # handle precedence for the people who don't use parentheses
            while ops and ops[-1] != "(" and precedence.get(ops[-1],0) >= precedence[token]:
                ops, output = __apply_op(ops, output)
            ops.append(token)
        else:
            # a simple condition string, e.g. "mass1 > 4"
            output.append(_parse_condition(token))

    # apply remaining ops
    while ops:
        if ops[-1] == "(":
            raise ValueError("Mismatched parentheses")
        ops, output = __apply_op(ops, output)

    if len(output) != 1:
        raise ValueError("Invalid expression")
    return output[0]


def select_phase(bpp, condition):
    '''Select the rows of the BPP array corresponding to a given phase.
    
    Parameters
    ----------
    bpp : pd.DataFrame
        DataFrame of the BPP array from COSMIC
    condition : str or Group or Condition, optional
        The condition by which to select the phase. This can be a pre-defined phase name (one of
        ["BBH_merger", "BNS_merger", "NSBH_merger", "BH_MS", "NS_MS", "WD_MS", "BH_GS", "NS_GS", "WD_GS"]), or
        a custom condition string (e.g. "mass_1 > 4 & mass_2 < 3") based on the column names of the BPP
        DataFrame. These conditions can have nested parentheses and use the operators &, |, and ~ for
        and, or, and not, respectively. The comparison operators allowed are ==, !=, <, <=, >, >=, "in",
        and "not in". See examples below for more details.    
    
    Returns
    -------
    out : pd.DataFrame or None
        DataFrame at the time of the selected phase

    Examples
    --------
    >>> # select all binary neutron star mergers
    >>> bns = select_phase(bpp, 'BNS_merger')

    >>> # select all binaries where either star is a main sequence star with teff > 5000 K
    >>> ms_hot = select_phase(bpp, '(teff_1 > 5000 & kstar_1 in [0,1]) | (teff_2 > 5000 & kstar_2 in [0,1])')

    >>> # select all binaries with eccentricity > 0.9 and separation < 10 Rsun
    >>> ecc_sep = select_phase(bpp, 'ecc > 0.9 & sep < 10')

    >>> # select all binaries that have undergone a supernova while unbound
    >>> sn_unbound = select_phase(bpp, 'evol_type in [15,16] and sep < 0')
    '''
    # handle custom conditions and groups
    if isinstance(condition, str):
        if condition in DEFAULT_PHASES:
            condition = DEFAULT_PHASES[condition]
        condition = _parse_conditions(condition)

    mask = _group_to_mask(bpp, condition)
    return bpp.loc[mask]
