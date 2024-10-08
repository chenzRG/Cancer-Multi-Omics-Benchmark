#!/bin/bash

DROPBOX_LINK="https://www.dropbox.com/scl/fi/g2djwgblekfepd62kgria/Main_Dataset.zip?rlkey=l5rl49qqnic6tpc9b3hyzjbmb&st=j4c7g5lz&dl=1"

# Destination directory for the download
DEST_DIR=".."

# Extract the file name from the Dropbox link
FILE_NAME=$(basename "$DROPBOX_LINK" | cut -d'?' -f1)

# Modify the download link (change dl=0 to dl=1)
DOWNLOAD_LINK="${DROPBOX_LINK/dl=0/dl=1}"

# Download the file to the specified directory
curl -L -o "${DEST_DIR}/${FILE_NAME}" "$DOWNLOAD_LINK"

echo "File has been downloaded to ${DEST_DIR}/${FILE_NAME}."

# Unzip the downloaded folder
unzip "${DEST_DIR}/${FILE_NAME}" -d "$DEST_DIR"

# Remove the ZIP file after extraction
rm "${DEST_DIR}/${FILE_NAME}"

echo "Folder has been downloaded and extracted to ${DEST_DIR}/${DEST_DIR}."