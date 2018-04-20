import argparse

from jetstream.utils import write_test_data


def main(args=None):
    parser = argparse.ArgumentParser(
        description="Create a table with test data.")
    parser.add_argument('path', help='Path to write out table')
    parser.add_argument(
        '--dialect',
        default='unix',
        choices=('unix', 'excel', 'excel-tab'),
        help='Dialect to use, see '
    )

    args = parser.parse_args(args)
    print(args)

    path = args.path
    dialect = args.dialect

    write_test_data(path, dialect=dialect)
    print('Wrote test data to {}'.format(path))


if __name__ == "__main__":
    main()
