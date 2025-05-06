input_list = [
    ['01', [
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16777734.0']],
    ['02', [
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16778265.0',
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16778265.1'
    ]],
    ['03', [
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16778612.0',
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16778612.1',
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16778612.1',
    ]],
    ['04', [
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16777915.0',
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16777915.1',
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16777915.2',
        '/home/b/b383366/sky/scc_at_isc25/nvsmi.log.16777915.3',
    ]],
]

def read_data_from_nvsmi_file(filename):
    """
    :param str filename:                Filename of the log file created with nvidia-smi

    :returns:                           List with two entries: \
                                        1) List with time data (x-axis) in [s] without offset, and \
                                        2) Dictionary with data read for GPUs 1-4
    :rtype:                             list[bool, list[float], dict]
    """
    import datetime

    ref_timestamp = None
    x_time = []
    y_gpu = {
        '0': [],
        '1': [],
        '2': [],
        '3': []
    }

    with open(filename, 'r') as infile:
        for single_line in infile:
            single_line = single_line.strip()
            new_stamp = False
            try:
                timestamp = datetime.datetime.strptime(single_line[:26], '%Y-%m-%dT%H:%M:%S,%f')
                new_stamp = True
            except Exception:
                pass

            if (new_stamp):
                if (ref_timestamp is None):
                    x_time = [0]
                    ref_timestamp = timestamp
                else:
                    delta_t = timestamp - ref_timestamp
                    x_time += [delta_t.seconds + 0.000001*delta_t.microseconds]

            elif any([single_line.startswith(x + ',') for x in y_gpu.keys()]):
                parts = single_line.split(',')
                y_gpu.update({parts[0]: y_gpu[parts[0]] + [[float(x.strip().split()[0]) for x in parts[1:]]]})

    return [x_time, y_gpu]


def process_files(filenamelist):
    data = {}
    for filename in filenamelist:
        x_time, y_gpu = read_data_from_nvsmi_file(filename)
        data.update({
            filename: {
                'x_time': x_time,
                'y_gpu': y_gpu
            }
        })

    return data
    

def average_runtime_and_energy(filenamelist):
    data = process_files(filenamelist=filenamelist)
    avg_energy = []
    runtimes = []
    for filename, filedata in data.items():
        duration = data[filename]['x_time'][-1]
        gpulist = sorted(list(filedata['y_gpu'].keys()))
        for gpu in gpulist:
            gpudata = filedata['y_gpu'][gpu]
            if (gpudata == []):
                continue

            tmp_avg = [sum([x[idx_value] for x in gpudata])/len(gpudata) \
                for idx_value in range(len(gpudata[0]))]
            avg_energy += [tmp_avg[0]*duration]
            runtimes += [duration]

    return [runtimes, avg_energy]


results = []
for nodesize, filenamelist in input_list:
    runtimes, avg_energy = average_runtime_and_energy(filenamelist)

    # Average energy used on all nodes
    print('Results for configuration ' + nodesize)
    #print('Energy over runtime: ' + str(avg_energy))
    print('(max.) time (s): ' + str(max(runtimes)))
    energy_used = sum(avg_energy)/3600
    print('Energy used (Wh): ' + str(energy_used))
    print('')
    results += [[int(nodesize), max(runtimes), energy_used]]

# results = [
#     [1, 389.107228, 60.11344661480575],
#     [3, 500.094681, 178.21741853835195],
#     [4, 659.684862, 288.51646378530705]
# ]


def plot_runtime_energy(results):
    from matplotlib import pyplot

    nodes = [x[0] for x in results]
    runtimes = [x[1] for x in results]
    energy = [x[2] for x in results]

    col1 = [0.6, 0.6, 0.8]
    col2 = [0.9, 0.5, 0.5]
    fig = pyplot.figure()
    ax1 = fig.gca()
    ax1.grid(True, which='both', linestyle=':')

    ax1.set_xlabel('Nodes used')
    ax1.plot(nodes, runtimes, color=col1)
    ax1.tick_params('y', colors=col1)
    ax1.set_ylabel('Runtime in s', color=col1)
    ax1.set_ylim(350, 700)

    ax2 = ax1.twinx()
    ax2.plot(nodes, energy, color=col2)
    ax2.tick_params('y', colors=col2)
    ax2.set_ylabel('Energy in Wh', color=col2)
    ax2.set_ylim(0, 350)

    #pyplot.show()
    pyplot.savefig("Energy.pdf")


def plot_speedup(results):
    from matplotlib import pyplot

    nodes = [x[0] for x in results]
    relative_runtimes = [x[1]/results[0][1] for x in results]
    relative_energy = [x[2]/results[0][2] for x in results]

    fig = pyplot.figure()
    ax = fig.gca()
    ax.grid(True, which='both', linestyle=':')

    ax.set_xlabel('Nodes used')
    ax.plot(nodes, [1/x for x in relative_runtimes])
    ax.tick_params('y')
    ax.set_ylabel('Speedup')

    #pyplot.show()
    pyplot.savefig("Speedup.pdf")


def plot_energy_speedup(results):
    from matplotlib import pyplot

    nodes = [x[0] for x in results]
    relative_runtimes = [x[1]/results[0][1] for x in results]
    relative_energy = [x[2]/results[0][2] for x in results]

    col1 = [0.6, 0.6, 0.8]
    col2 = [0.9, 0.5, 0.5]
    fig = pyplot.figure()
    ax = fig.gca()
    ax.grid(True, which='both', linestyle=':')

    ax.set_xlabel('Nodes used')
    ax.plot(nodes, [x*y for x, y in zip(relative_energy, relative_runtimes)])
    ax.tick_params('y')
    ax.set_ylabel('Relative Energy/Speedup')

    #pyplot.show()
    pyplot.savefig("RelativeEnergy.pdf")


plot_runtime_energy(results)
plot_speedup(results)
plot_energy_speedup(results)
