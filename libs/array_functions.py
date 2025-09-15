import numpy as np
import random

def find_nearest_index(array,value):

    '''
    inputs:
        array (1D array or slice of higher dimensional array) -- the array to search for the value
        value (float) -- the value to search for
    outputs:
        idx (int) -- the index of array element whose value is closest to the search value
    '''

    idx = (np.abs(array - value)).argmin()
    return idx

def get_indices_within(array,value,threshold):

    '''
    inputs:
        array (1D array or slice of higher dimensional array) -- the array to search for the value
        value (float) -- the value to search for
        threshold (float) -- acceptable deviation from the search value
    outputs:
        arr (1D array of ints) -- array containing the indices of array elements whose values are less than the threshold away from the search value.
    '''

    arr = []
    for index,x in np.ndenumerate(array):
        if np.abs(x - value) <= threshold:
            arr.append(index)
    return np.array(arr)

def parabolic_max(xs,ys):

    '''
    inputs:
        xs (1D array containing 3 or fewer elements) -- x coordinates of points you want to interpolate
        ys (1D array containing 3 or fewer elements) -- y coordinates of points you want to interpolate
    outputs:
        max_tuple (2-tuple of floats) -- if 3 points are provided, contains the coords of the maximum (or minimum) specified by these points.
            if < 3 points are provided, averages of the coordinates of each point are returned.
    '''

    if xs[0] != xs[1] and xs[0] != xs[2] and xs[1] != xs[2]:
        x_mat = np.array([[xs[0]**2,xs[0],1],
                        [xs[1]**2,xs[1],1],
                        [xs[2]**2,xs[2],1]])
        x_mat_inv = np.linalg.inv(x_mat)
        y_vec = np.array([[ys[0]],
                        [ys[1]],
                        [ys[2]]])
        sol = np.matmul(x_mat_inv,y_vec)
        a = sol[0]
        b = sol[1]
        c = sol[2]
        return(-b/(2*a),c - b**2/(4*a))
    else:
        return((1/3)*(xs[0]+xs[1]+xs[2]),(1/3)*(ys[0]+ys[1]+ys[2]))

def generate_random_list(total_sum, num_elements):
    """
    Generates a list of random positive integers with a specified sum and number of elements.
    
    Args:
        total_sum (int): The target sum of the list elements.
        num_elements (int): The number of elements in the list.

    Returns:
        list: A list of positive integers whose sum equals total_sum.
    """
    if num_elements <= 0:
        return []
    if num_elements == 1:
        return [total_sum]

    # Generate num_elements - 1 random "breaks"
    breaks = sorted([random.uniform(1, total_sum - 1) for _ in range(num_elements - 1)])

    # Add 0 at the start and total_sum at the end
    breaks = [0] + breaks + [total_sum]

    # Calculate the differences between consecutive breaks
    result = [breaks[i] - breaks[i-1] for i in range(1, len(breaks))]
    
    return result

