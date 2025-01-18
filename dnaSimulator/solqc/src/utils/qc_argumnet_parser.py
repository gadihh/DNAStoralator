from argparse import ArgumentParser
import time
from datetime import datetime

def parse_arguments(arguments):
    parser = ArgumentParser(description="Possible arguments for the Quality-Control tool.")

    parser.add_argument('-d', '-design', dest="design", type=str)
    parser.add_argument('-r', '-reads', dest="reads", type=str)
    parser.add_argument('-c', '-config', dest="config", type=str)
    parser.add_argument('-m', '-aligners', dest="aligners", type=str, nargs='+', default=['BarcodeAligner', 'EditDistanceAligner'])
    parser.add_argument('-a', '-analyzers', dest='analyzers', type=str, nargs='+', default=['FrequencyAnalyzer'])
    parser.add_argument('-time', '-time', dest='time', type=str, nargs='?',
                        default=datetime.utcfromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S'))

    parser.add_argument('-id', dest='id', type=str, nargs='?', default="")

    parser.add_argument('--override', dest='override', action='store_true')
    parser.add_argument('--no-override', dest='override', action='store_false')
    parser.set_defaults(override=True)

    parser.add_argument('--edit', dest='edit', action='store_true')
    parser.add_argument('--no-edit', dest='edit', action='store_false')
    parser.set_defaults(edit=False)

    args = parser.parse_args(arguments)

    return args
