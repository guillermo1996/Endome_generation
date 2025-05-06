#!/bin/bash

echo "This is stdout output"
echo "This is stderr output" >&2
ls /nonexistent_directory  # This will generate an error message (stderr)
echo "Another stdout message"
