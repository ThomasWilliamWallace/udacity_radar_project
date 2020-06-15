IMPLEMENTATION STEPS FOR THE 2D CFAR PROCESS:
Start with a 2d matrix, the result of a 2d FFT.
For each cell
    If the cell is far enough from the edges that it has a full set of training and guard cells then
        Sum the db2pow values of the surrounding training cells, and take the average.
        Then convert this average back into dB with pow2db.
        This is the noise level of this cell
        If this cell value is greater than the noise level + a threshold, then
            Assign the cell a CFAR response of 1
        Else
            Assign the cell a CFAR response of 0
    Else
        Assign the cell a CFAR response of 0


For each cell, calculate the average of it's surrounding training cells, while excluding the guard cells, and the cell-under-test.

SELECTION OF TRAINING, GUARD CELLS AND OFFSET:
I chose 4 guard cells. When the cell-under-test is centered on the signal, this excludes the bulk of the signal, so that we are just measuring the surrounding noise.

I chose 3 training cells. This is the smallest number which produces a smooth averaged noise.
It's not any bigger as we still want it to reflect the local noise around the cell-under-test.

The offset of 2 was chosen because in this exercise it was the smallest offset which still comfortably filtered out the noise.
I chose to minimise this offset because a larger offset loses sensivity.

STEPS TAKEN TO SUPPRESS THE NON-THRESHOLDED CELLS AT THE EDGES:
I only calculated the noise, and CFAR response, for cells-under-test which had a full set of training and guard cells.
Any cell too close to the edge to have these, was instead assigned a CFAR response of zero.
