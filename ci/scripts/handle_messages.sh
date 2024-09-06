#!/usr/bin/env bash

# Check if tput is available
if command -v tput &> /dev/null; then
    # Define ANSI escape codes for colors and styles using tput
    bold=$(tput bold)
    red=$(tput setaf 1)
    green=$(tput setaf 2)
    yellow=$(tput setaf 3)
    reset=$(tput sgr0)
else
    # Disable warnings in non-interactive shells (e.g., CI environments)
    if [[ $- == *i* ]]; then
        echo "[WARNING] tput not found, colors and styles won't be applied."
    fi
    bold=""
    red=""
    green=""
    yellow=""
    reset=""
fi

# Function to print messages with specified color
# Arguments:
#   $1: Color code
#   $2: Message type (INFO, WARNING, ERROR)
#   $3: The actual message to be printed
print_message() {
    local color=$1
    local type=$2
    local message=$3
    echo -e "${bold}${color}[${type}]${reset} ${message}"
}

# Function to print info message
# Arguments:
#   $1: The info message to be printed
info() {
    local message=$1
    print_message "$green" "INFO" "$message"
}

# Function to print warning message
# Arguments:
#   $1: The warning message to be printed
warning() {
    local message=$1
    print_message "$yellow" "WARNING" "$message"
}

# Function to print error message
# Arguments:
#   $1: The error message to be printed
# Outputs the message to stderr
error() {
    local message=$1
    print_message "$red" "ERROR" "$message" >&2
}

# Function to handle errors
# Arguments:
#   $1: The error message to be printed before exiting
error_exit() {
    local message=$1
    error "$message"
    exit 1
}

# Handle any unexpected error in the script and exit gracefully
trap 'error_exit "An unexpected error occurred."' ERR
