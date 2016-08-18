"""
Python wrapper for interacting with hybrid_ss_min.
"""

import os
import pty
from subprocess import Popen, PIPE


class HybridSSMinWrapper(object):
    """An object that opens a hybrid_ss_min process and allows repeatedly
    sending it data without having to open a new process each time.
    """

    # CMD = HYBRID_SS_MIN_BINARY + ' --stream'
    CMD = 'hybrid-ss-min --stream'

    def __init__(self, extra_options=''):
        """Opens the process.
        """
        self.cmd = self.CMD + ' ' + extra_options
        self.master, self.slave = pty.openpty()
        self.process = Popen(
                self.cmd, shell=True, stdin=PIPE, stdout=self.slave,
                close_fds=True)
        self.stdout = os.fdopen(self.master)

    def compute_delta_g(self, sequence_string):
        """Computes the delta_g value for the given sequence_string.

        Returns:
            A float representing the free energy.
        """
        self.process.stdin.write(sequence_string + ';')
        return float(self.stdout.readline().strip())


if __name__ == '__main__':
    wrapper_obj = HybridSSMinWrapper()
    print wrapper_obj.compute_delta_g('GGGGGGCCCCC')
    print wrapper_obj.compute_delta_g('AAAAAAAAATTTTTTTTTTT')
    print wrapper_obj.compute_delta_g('AAAAAAAAAAAAAAAAAAAAAAAAAAAA')
