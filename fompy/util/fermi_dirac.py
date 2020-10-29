"""
This module contains functions for calculating the Fermi-Dirac integral.
The original code was taken from https://github.com/scott-maddox/fdint and translated from cython to python.
Unfortunately, the changes made the code about 35 times slower.
"""

# Copyright (c) 2015, Scott J Maddox.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
# disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
# following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
# products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from math import exp, sqrt


def fd1(phi):
    """
    Calculate the Fermi-Dirac integral.

    Parameters
    ----------
    phi : float

    Returns
    -------
    float
        The Fermi-Dirac integral of `phi`.
    """
    if phi < -2e0:
        return _fd1h_lt_m2(phi)
    if phi < 0e0:
        return _fd1h_m2_to_0(phi)
    if phi < 2e0:
        return _fd1h_0_to_2(phi)
    if phi < 5e0:
        return _fd1h_2_to_5(phi)
    if phi < 10e0:
        return _fd1h_5_to_10(phi)
    if phi < 20e0:
        return _fd1h_10_to_20(phi)
    if phi < 40e0:
        return _fd1h_20_to_40(phi)
    return _fd1h_gt_40(phi)


# @formatter:off


def _fd1h_lt_m2(phi):
    exp_phi = exp(phi)
    t = exp_phi * 7.38905609893065023e0
    return exp_phi * (0.886226925452758014e0
        - exp_phi * (19894.4553386951666e0
        + t * (4509.64329955948557e0
        + t * (303.461789035142376e0
        + t * (5.7574879114754736e0
        + t * 0.00275088986849762610e0
        )))) / (63493.915041308052e0
        + t * (19070.1178243603945e0
        + t * (1962.19362141235102e0
        + t * (79.250704958640158e0
        + t)))))


def _fd1h_m2_to_0(phi):
    s = -0.5e0 * phi
    t = 1e0-s
    return (149.462587768865243e0
        + t * (22.8125889885050154e0
        + t * (-0.629256395534285422e0
        + t * (9.08120441515995244e0
        + t * (3.35357478401835299e0
        + t * (-0.473677696915555805e0
        + t * (-0.467190913556185953e0
        + t * (-0.0880610317272330793e0
        -t * 0.00262208080491572673e0
        )))))))) / (269.94660938022644e0
        + s * (343.6419926336247e0
        + s * (323.9049470901941e0
        + s * (218.89170769294024e0
        + s * (102.31331350098315e0
        + s * (36.319337289702664e0
        + s * (8.3317401231389461e0
        + s)))))))


def _fd1h_0_to_2(phi):
    t = 0.5e0 * phi
    return (71652.717119215557e0
        + t * (134954.734070223743e0
        + t * (153693.833350315645e0
        + t * (123247.280745703400e0
        + t * (72886.293647930726e0
        + t * (32081.2499422362952e0
        + t * (10210.9967337762918e0
        + t * (2152.71110381320778e0
        + t * 232.906588165205042e0
        )))))))) / (105667.839854298798e0
        + t * (31946.0752989314444e0
        + t * (71158.788776422211e0
        + t * (15650.8990138187414e0
        + t * (13521.8033657783433e0
        + t * (1646.98258283527892e0
        + t * (618.90691969249409e0
        + t * (-3.36319591755394735e0
        + t))))))))


def _fd1h_2_to_5(phi):
    t = 0.3333333333333333333e0 * (phi-2e0)
    return (23744.8706993314289e0
        + t * (68257.8589855623002e0
        + t * (89327.4467683334597e0
        + t * (62766.3415600442563e0
        + t * (20093.6622609901994e0
        + t * (-2213.89084119777949e0
        + t * (-3901.66057267577389e0
        -t * 948.642895944858861e0
        ))))))) / (9488.61972919565851e0
        + t * (12514.8125526953073e0
        + t * (9903.44088207450946e0
        + t * (2138.15420910334305e0
        + t * (-528.394863730838233e0
        + t * (-661.033633995449691e0
        + t * (-51.4481470250962337e0
        + t)))))))


def _fd1h_5_to_10(phi):
    t = 0.2e0 * phi-1e0
    return (311337.452661582536e0
        + t * (1.11267074416648198e6
        + t * (1.75638628895671735e6
        + t * (1.59630855803772449e6
        + t * (910818.935456183774e0
        + t * (326492.733550701245e0
        + t * (65507.2624972852908e0
        + t * 4809.45649527286889e0
        ))))))) / (39721.6641625089685e0
        + t * (86424.7529107662431e0
        + t * (88163.7255252151780e0
        + t * (50615.7363511157353e0
        + t * (17334.9774805008209e0
        + t * (2712.13170809042550e0
        + t * (82.2205828354629102e0
        - t))))))) * 0.999999999999999877e0


def _fd1h_10_to_20(phi):
    t = 0.1e0 * phi-1e0
    return (7.26870063003059784e6
        + t * (2.79049734854776025e7
        + t * (4.42791767759742390e7
        + t * (3.63735017512363365e7
        + t * (1.55766342463679795e7
        + t * (2.97469357085299505e6
        + t * 154516.447031598403e0
        )))))) / (340542.544360209743e0
        + t * (805021.468647620047e0
        + t * (759088.235455002605e0
        + t * (304686.671371640343e0
        + t * (39289.4061400542309e0
        + t * (582.426138126398363e0
        + t * (11.2728194581586028e0
        - t)))))))


def _fd1h_20_to_40(phi):
    t = 0.05e0 * phi-1e0
    return (4.81449797541963104e6
        + t * (1.85162850713127602e7
        + t * (2.77630967522574435e7
        + t * (2.03275937688070624e7
        + t * (7.41578871589369361e6
        + t * (1.21193113596189034e6
        + t * 63211.9545144644852e0
        )))))) / (80492.7765975237449e0
        + t * (189328.678152654840e0
        + t * (151155.890651482570e0
        + t * (48146.3242253837259e0
        + t * (5407.08878394180588e0
        + t * (112.195044410775577e0
        - t))))))


def _fd1h_gt_40(phi):
    factor = 2e0/3e0   # = 1/(k + 1)
    w = 1e0/(phi * phi)
    s = 1e0-1600e0 * w
    return phi * sqrt(phi) * factor * (1e0 + w
        * (8109.79390744477921e0
        + s * (342.069867454704106e0
        + s * 1.07141702293504595e0))
        / (6569.98472532829094e0
        + s * (280.706465851683809e0
        + s)))

# @formatter:on
