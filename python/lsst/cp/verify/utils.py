# This file is part of cp_verify.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
__all__ = ['mergeStatDict']


def mergeStatDict(statDictA, statDictB):
    """Merge one dictionary of statistics with another.

    Parameters
    ----------
    statDictA, statDictB : `dict` [`str`, `dict`]
        Input dictionaries of statistics.

    Returns
    -------
    output : `dict` [`str`, `dict`]
        A dictionary containing the union of the internal
        dictionaries of ``statDictA`` and ``statDictB``, indexed
        by the amplifier names.

    Raises
    ------
    RuntimeError :
        Raised if the merge would overwrite a key/value pair, or if
        the amplifier names do not match between the two inputs.
    """
    if set(statDictA.keys()) != set(statDictB.keys()):
        raise RuntimeError(f"Mismatch in amplifier names: {statDictA.keys()} {statDictB.keys()}")

    for amp in statDictA.keys():
        if not statDictA[amp].keys().isdisjoint(statDictB[amp].keys()):
            raise RuntimeError(f"Duplicate keys passed: {statDictA[amp].keys()} {statDictB[amp].keys()}")
        statDictA[amp].update(statDictB[amp])
    return statDictA
