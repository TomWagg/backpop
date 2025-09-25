from argparse import ArgumentParser
from .backpop import BackPop

def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--ini_file', help='Path to INI file', type=str, default="params.ini")
    args = parser.parse_args()

    bp = BackPop(config_file=args.ini_file)
    bp.run_sampler()
    bp.save_output()

if __name__ == "__main__":
    main()
