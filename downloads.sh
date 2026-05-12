#!/bin/bash 
# This script is used to download the data from the CSV files and save them in the data directory.
# It uses the wget command to download the files from the URLs specified in the urls.txt file
# The downloaded files are saved in the data directory with the same name as the URL
# The script also checks if the data directory exists, if not it creates it before downloading the files
# Create data directory if it doesn't exist
mkdir -p data
# Read URLs from urls.txt and download each file
while IFS= read -r url; do

    # Extract filename from URL
    filename=$(basename "$url")
    # Download the file using wget
    wget -O "data/$filename" "$url"
done < urls.txt

