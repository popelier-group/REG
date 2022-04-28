"""
setup_rdp.py v0.1
F. Falcioni, P. L. A. Popelier

Script with the function to run the Ramer-Douglas-Peucker algorithm to potential energy surfaces.

Check for updates at github.com/FabioFalcioni

Please, report bugs and issues to fabio.falcioni@manchester.ac.uk
coded by F.Falcioni
"""
import reg
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser


def get_energies_from_file(file, file_type):
    '''
    Function to parse .txt files containing energies at each step of a PES
    '''
    energies = []
    if file_type == "txt":
        with open(file, 'r') as f:
            lines = f.readlines()[4:]
            for line in lines:
                for value in (line.split()[1:]):
                    energies.append(float(value))
    elif file_type == "grep":
        with open(file, 'r') as f:
            for line in f.readlines():
                energies.append(float(line.split(' ')[-1]))
    return energies


def find_point_between_two_points(start, end, value):
    '''
    Function to find the y value of a point between two other points of a function given its x value (i.e. interpolation)
    '''
    slope = slope_intercept(start[0], start[1], end[0], end[1])
    intercept = start[1] - slope * start[0]
    y = slope * value + intercept
    return y


def slope_intercept(x1, y1, x2, y2):
    '''
    Function to find the slope and intercept of a line between two points
    '''
    a = (y2 - y1) / (x2 - x1)
    return a


def deviation(point, start, end):
    '''
    Function to calculate the perpendicular distance of a point from the line defined by other two points.
    '''
    if np.all(np.equal(start, end)):
        return np.linalg.norm(point - start)

    return np.divide(
        np.abs(np.linalg.norm(np.cross(end - start, start - point))),
        np.linalg.norm(end - start))


def rdp_rec(M, epsilon, pdist=deviation):
    '''
    Function to run the RDP algorithm recursively
    '''
    max_dist = 0.0
    index = -1

    for i in range(1, M.shape[0]):
        d = pdist(M[i], M[0], M[-1])

        if d > max_dist:
            index = i
            max_dist = d

    if max_dist > epsilon:
        r1 = rdp_rec(M[:index + 1], epsilon, pdist)
        r2 = rdp_rec(M[index:], epsilon, pdist)

        return np.vstack((r1[:-1], r2))
    else:
        return np.vstack((M[0], M[-1]))


def rdp(M, epsilon=0, pdist=deviation):
    if "numpy" in str(type(M)):
        return rdp_rec(M, epsilon, pdist)

    return rdp_rec(np.array(M), epsilon, pdist).tolist()


def find_maximum_deviation(M, pdist=deviation):
    '''
    Function to calculate the maximum perpendicular distance from the line defined between the extreme points of a segment
    '''
    max_dist = 0.0
    index = -1

    for i in range(1, M.shape[0]):
        d = pdist(M[i], M[0], M[-1])

        if d > max_dist:
            index = i
            max_dist = d
    return max_dist


def rmse(true_value, predicted_value):
    '''
    Function to calculate the Root-Mean-Squared-Error between true and predicted values
    '''
    if len(true_value) != len(predicted_value):  ## Checks if X and Y have the same size
        raise ValueError("Arrays must have the same size")
    error = [true_value[i] - predicted_value[i] for i in range(len(true_value))]
    temp = [a ** 2 for a in error]
    RMSE = (sum(temp) / len(temp)) ** 0.5
    return RMSE


def minimum_points_segment_RDP(energy, cc, epsilon):
    '''
    Function to obtain a segment (or function) with the minimum amount of points given a specific value of epsilon (i.e. deviation)
    '''
    vector = []
    x = []
    y = []

    for i in range(0, len(energy)):
        coordinate = (cc[i], energy[i])
        vector.append(coordinate)
    new_points = rdp(vector, epsilon=epsilon)
    for value in new_points:
        new_x, new_y = value
        x.append(new_x)
        y.append(new_y)
    new_interpolated_points = []
    for j in range(1, len(new_points)):
        for i in range(0, len(cc)):
            if cc[i] > new_points[j - 1][0] and cc[i] < new_points[j][0]:
                interpolate_point = find_point_between_two_points(new_points[j - 1], new_points[j], cc[i])
                new_interpolated_points.append(interpolate_point)
    y_reference = y + new_interpolated_points
    if energy[0] > energy[-1]:
        y_reference.sort(reverse=True)
    else:
        y_reference.sort(reverse=False)
    RMSE = rmse(energy, y_reference)
    return y, x, RMSE


def cross_validation_RDP(energy, cc, rmse_confidence=1, step_size=0.01):
    '''
    Function to run the RDP algorithm X times (depending on step_size) at different values of epsilon and obtain
    obtain a function with lowest number of points at that (or below) RMSE of confidence.
    '''
    new_vector = []
    [new_vector.append((cc[i], energy[i])) for i in range(0, len(energy))]
    max_epsilon = find_maximum_deviation(np.array(new_vector), pdist=deviation)
    epsilon = np.arange(0, max_epsilon, step_size).tolist()
    y_values_minimum = []
    x_values_minimum = []
    rmse_min = []
    for eps in epsilon:
        y_values, x_values, rmse = minimum_points_segment_RDP(energy, cc, eps)
        if rmse <= rmse_confidence:
            y_values_minimum.append(y_values)
            x_values_minimum.append(x_values)
            rmse_min.append(rmse)

    return min(x_values_minimum, key=len), min(y_values_minimum, key=len), rmse_min[
        x_values_minimum.index(min(x_values_minimum, key=len))]


def main():
    #Parsing USER INPUT
    global segments, cc_new
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-f", "--file", action='store', type='string', dest='file', help="input file with PES energies")
    parser.add_option("-r", "--rmse", action='store', type='float', dest='rmse_confidence',
                      help="Root Mean Squared error of confidence")
    parser.add_option("-c", "--critical_points", action='store', type='string', dest='critical_points',
                      help="'AUTO' for automatic search / List of integer values (e.g. cp1,cp2,..cpn) for user "
                           "defined points / Nothing for single segment PES")
    (option, args) = parser.parse_args()

    print(
        "RDP setup: searching for a new polyline with RMSE of confidence {} kJ/mol ...".format(option.rmse_confidence))
    wfn_energies = get_energies_from_file(option.file, file_type='txt')
    wfn_energies = np.array(wfn_energies) * 2625.5 #Converting in kJ/mol from a.u.
    cc = np.array([i for i in range(1, len(wfn_energies) + 1)])

    # Search for critical points in the function
    if option.critical_points == 'AUTO':
        tp = reg.find_critical(wfn_energies, cc, use_inflex=False, min_points=5)
        segments = reg.split_segm(wfn_energies - sum(wfn_energies) / len(wfn_energies), tp)
        cc_new = reg.split_segm(cc, tp)
    elif isinstance(option.critical_points, str):
        tp = [(int(value)) for value in option.critical_points.split(',')]
        print(tp)
        segments = reg.split_segm(wfn_energies - sum(wfn_energies) / len(wfn_energies), tp)
        cc_new = reg.split_segm(cc, tp)
    elif not option.critical_points:
        segments = [wfn_energies - sum(wfn_energies) / len(wfn_energies)]
        cc_new = [cc]

    #Creating a figure with the new polyline given the USER input rmse of confidence
    plt.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=( 9,5))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    for i in range(0, len(segments)):
        x, y, RMSE = cross_validation_RDP(segments[i], cc_new[i], rmse_confidence=option.rmse_confidence)
        print("RMSE Segment {} = {}".format(i + 1, RMSE))
        print("Points of the new polyline for Segment {} = {}".format(i + 1, x))
        plt.plot(x, y, marker='o', markersize=10)
        plt.plot(x[0], y[0], marker='o', markersize=10, c='red')
        plt.plot(x[-1], y[-1], marker='o', markersize=10, c='red')
        textstr = 'RMSE = {}'.format(round(RMSE, 2))
        plt.text(sum(x) / len(x), sum(y) / len(y), textstr, fontsize=12, verticalalignment='top', bbox=props)
    plt.plot(cc, wfn_energies - sum(wfn_energies) / len(wfn_energies), c='#4d4d4d', marker='o', markersize=2)
    plt.xlabel(r'Control Coordinate (N step)')
    plt.ylabel(r'Relative Energy (kJ $\mathrm{mol^{-1}}$)')
    plt.show()
    fig.savefig('RDP_out.png', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    main()
