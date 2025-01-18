from PyQt5.QtCore import *
import ast
from index_clustering import Clustering
from hash_based_clustering import *
from simulator import Simulator
from filepath import *
from result import *


# abstract base class for workers:
class Worker(QObject):
    """
    According to the documentation about threading and parallelism, it is preferred to wrap a thread by a QObject class.
    Inside this class the QThread will be created and using a method moveToThread it will be incorporated into a wrapper.
    Two workers are subclassed: ErrorSimWorker and ClusteringWorker, they all have a simulator object which performs
    the main loop and inside this loop a progress update is being sent as a signal as well as a termination flag checked
    for the need to cancel the simulation. When a stop button pressed for a specific worker, this flag will be set on.
    """
    update_progress = pyqtSignal(int)
    work_done = pyqtSignal()
    terminated = pyqtSignal()

    def __init__(self, path):
        super().__init__()
        self.input_path = path
        self.simulator = None
        self.thread = QThread()
        self.moveToThread(self.thread)
        self.thread.started.connect(self.work)
        self.work_done.connect(self.thread.quit)

        self.files = []

    def report_func(self, total_lines, curr_line, to_terminate):
        percent = int(curr_line * 100 // total_lines)
        self.update_progress.emit(percent)
        if to_terminate:
            raise ThreadTerminated

    def start(self):
        self.thread.start()

    def kill(self):
        # inside simulator's main loop there is a constant check for this field. if set - termination will happen.
        self.simulator.to_terminate = True
        self.thread.quit()

    def work(self):
        raise NotImplementedError()


class SimulateErrorsWorker(Worker):

    def __init__(self, general_errors, per_base_errors, input_dna_path, stutter_chosen, dist_info, error_output_path,
                 shuffled_output_path):
        # error_output_path will be stored in input_path for generalization of Worker super class
        super(SimulateErrorsWorker, self).__init__(error_output_path)
        self.general_errors = general_errors
        self.per_base_errors = per_base_errors
        self.inputDNAPath = input_dna_path
        self.stutter_chosen = stutter_chosen
        self.dist_info = dist_info
        self.shuffled_output_path = shuffled_output_path
        self.temp_file = remove_extension(self.input_path) + '_temp.txt'

        # parse value to list if type is vector:
        if self.dist_info is not None and self.dist_info['type'] == 'vector':
            self.dist_info['value'] = ast.literal_eval(self.dist_info['value'])
        self.simulator = Simulator(self.general_errors, self.per_base_errors, self.inputDNAPath, self.stutter_chosen,
                                   self.dist_info, self.input_path, self.shuffled_output_path, self.temp_file)

    def work(self):
        try:
            self.simulator.simulate_errors(self.report_func)
        except ThreadTerminated:
            # cancel button was pressed:
            self.terminated.emit()
            return
        else:
            # all work was done without interruptions:
            self.work_done.emit()


class ClusteringWorker(Worker):

    def __init__(self, cluster_tech, cluster_index, workdir, shuffled_path, orig_path, output_path):
        super(ClusteringWorker, self).__init__(orig_path)
        self.cluster_tech = cluster_tech
        self.cluster_index = cluster_index
        self.workdir = workdir
        self.input_shuffled_path = shuffled_path
        self.output_path = output_path

        self.simulator = Clustering(self.cluster_index, self.cluster_tech, self.workdir,
                                    self.input_path, self.input_shuffled_path, self.output_path)
        self.terminated.connect(self.simulator.cleanup)

    def work(self):
        start_time = time.time()
        self.simulator.create_all_keys()
        self.simulator.fill_dict_from_shuffled()
        try:
            self.simulator.full_edit_distance(self.report_func)
        except ThreadTerminated:
            self.terminated.emit()
            return
        else:
            self.simulator.create_orig_dict()
            [num_of_errors, num_of_false_negative] = self.simulator.compare_orig_with_clustering()

            clustering_result = ClusteringResult(self.input_path, self.cluster_index, self.cluster_tech,
                                                 round((time.time() - start_time), 2),
                                                 num_of_errors, num_of_false_negative,
                                                 str(self.simulator.total_amount_of_thrown))
            # call add to json:
            keys = dict(name=clustering_result.name, tech=clustering_result.tech, index=clustering_result.index)
            update_json_file(Clustering.res_file(), clustering_result.__dict__, keys)
            # signal the app that worker is done and data is ready.
            self.work_done.emit()

class HashBasedClusteringWorker(Worker):
    def __init__(self, cluster_tech, cluster_index, number_of_steps, windows_size, similarity_threshold, chosen_clustering_algo, clustering_result_dir, clustering_input, clustering_output):
        super(HashBasedClusteringWorker, self).__init__(clustering_input)
        self.cluster_index = cluster_index
        self.cluster_tech = cluster_tech
        self.chosen_clustering_algo = chosen_clustering_algo
        self.number_of_steps = number_of_steps
        self.windows_size = windows_size
        self.similarity_threshold = similarity_threshold
        self.results_f = clustering_result_dir + '_final_results.txt'
        self.clustering_input = clustering_input
        self.clustering_output = clustering_output
        self.simulator =  HashBasedCluster(self.cluster_tech)

    def work(self):
        start_time = time.time()
        self.simulator.hash_based_cluster(self.number_of_steps,  self.windows_size,
                                    self.similarity_threshold, self.report_func, self.clustering_input,self.clustering_output)

        total_clusters_count =  self.simulator.get_number_of_clusters()
        with open(self.results_f, 'w') as results_file:
            print("Run time: %s seconds" % round((time.time() - start_time),2), file=results_file)
            print('Number of clusters generated: ' + str(total_clusters_count), file=results_file)
        self.work_done.emit()


class ReconProcess(QProcess):
    def __init__(self, input_path, algo, parent=None):
        super().__init__(parent)
        self.input_path = input_path
        self.algo = algo


class ThreadTerminated(Exception):
    pass
