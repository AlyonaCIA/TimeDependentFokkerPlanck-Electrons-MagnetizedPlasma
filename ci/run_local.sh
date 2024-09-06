#!/usr/bin/env bash
## Parse arguments
set -e

POSITIONAL=()
# Set path to CI Scripts as environmental variable
export CI_SCRIPTS_PATH="ci/scripts"
export CI_FOLDER="ci/scripts"
# Change all shell scripts in the ci/ folder to have permission to execute
chmod u+x ${CI_FOLDER}/*.sh
# Set formatting and tests run to False initially
run_format=false
run_tests=false

function print_usage {
cat <<- EOF
    usage: run_local [-f | -t | -i | -h] -- Runs code quality and formatting checks

    -f,--run-format: Forces autoformatting in-place on Fortran and Python code.

    -t,--test,--tests: Runs format and code quality checks, along with unit tests.

    -i,--install: Install required packages for the environment.

    -h,--help: Get help on how to use the local code quality and test script.

EOF
}

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --run-format|-f)
        run_format=true
        shift
    ;;
    --tests|--test|-t)
        run_tests=true
        shift
    ;;
    --install|-i)
        echo "Installing requirements"
        # Install Python dependencies
        if [ -f "requirements.txt" ]; then
            pip install -r requirements.txt
        else
            echo "requirements.txt not found. Skipping Python dependencies installation."
        fi
        
        if [ -f "${CI_SCRIPTS_PATH}/requirements/validation_requirements.txt" ]; then
            pip install -r ${CI_SCRIPTS_PATH}/requirements/validation_requirements.txt
        else
            echo "Validation requirements file not found. Skipping validation dependencies installation."
        fi

        # Install Fortran dependencies (using apt for gfortran as an example)
        if command -v apt &> /dev/null; then
            sudo apt-get update
            sudo apt-get install -y gfortran  # Adjust this based on your requirements
        fi

        echo "Installation succeeded"
        exit 0
    ;;
    --help|-h)
        print_usage
        exit 0
    ;;
    *)
        echo "Unrecognized argument $key. Try ${CI_SCRIPTS_PATH}/run_local.sh --help"
        exit 2
    ;;
esac
done

if $run_format; then
    echo "Applying code formatting"
    # Apply Python formatting
    ${CI_SCRIPTS_PATH}/apply_formatting.sh

    # Apply Fortran formatting using fprettify
    if command -v fprettify &> /dev/null; then
        find . -name "*.f90" -o -name "*.f" | xargs fprettify --silent
    else
        echo "fprettify not found. Please ensure it is installed."
        exit 1
    fi
fi

# Run formatting validation and linting for both Python and Fortran
echo "Validating code formatting and running linting"
${CI_SCRIPTS_PATH}/requirements/validation_requirements.txt
${CI_SCRIPTS_PATH}/run_linting.sh

if $run_tests; then
    echo "Running tests"
    # Assuming pFUnit is used for unit testing in Fortran
    if command -v pfunit &> /dev/null; then
        mpirun -np 4 ./test_suite || { echo "Fortran tests failed."; exit 1; }
    else
        echo "pFUnit not found. Please ensure it is installed and available in your PATH."
        exit 1
    fi

    # Run Python tests if they exist
    if [ -f "tests/test_suite.py" ]; then
        pytest || { echo "Python tests failed."; exit 1; }
    else
        echo "No Python test suite found."
    fi
fi

echo "Local checks passed :) Might still deviate from GitHub Actions checks :/"
