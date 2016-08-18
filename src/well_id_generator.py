# Copyright (C) 2016 The President and Fellows of Harvard College
#
# This file may require additional software or modules to be installed to run
# run properly.
#
# The Genome Recoder software is available under an internal non-commercial
# research and academic use license.  Questions about this software or the
# licensing thereof can be addressed to  Office of Technology Development,
# Harvard University, email: otd@harvard.edu.
#
# @author Gleb Kuznetsov (kuznetsov@g.harvard.edu)

"""
Generator object for standard 96-well plate.
"""


class WellIdGenerator(object):
    """Generates 96-plate well ids from A1 ... A12, ..., H1 ... H12

    Also returns plate number if requested.
    """

    LETTER_TRANSITION_TABLE = {
            'A': 'B',
            'B': 'C',
            'C': 'D',
            'D': 'E',
            'E': 'F',
            'F': 'G',
            'G': 'H',
            'H': 'A',
    }


    def __init__(self, start_well_letter='A', start_well_number=1,
            start_plate=1, include_plate=False):
        assert start_well_letter in self.LETTER_TRANSITION_TABLE
        assert start_well_number > 0 and start_well_number <= 12
        self.letter = start_well_letter
        self.number = start_well_number

        # Track plate, but only show if included.
        self.plate = start_plate
        self.include_plate = include_plate

    def __iter__(self):
        return self

    def next(self):
        # Save current state for next return before updating.
        current_plate = self.plate
        current_id = self.letter + "%02d" % (self.number,)

        # Bump the state.
        if self.number == 12:
            self.number = 1
            self.letter = self.LETTER_TRANSITION_TABLE[self.letter]
            # If we are back to A, bump the plate number.
            if self.letter == 'A':
                self.plate += 1
        else:
            self.number += 1

        # Return current state.
        if self.include_plate:
            return (current_plate, current_id)
        return current_id
