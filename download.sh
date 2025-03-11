#!/bin/bash

# Google Drive file ID
FILE_ID="1JF4aHuRS7_EXgjN8VySej0wNfoXn9nnQ"

# Destination directory for the download
DEST_DIR=".."

# 文件名
FILE_NAME="Main_Dataset.zip"

# check gdown installation
if ! command -v gdown &> /dev/null
then
    echo "gdown not found, installing..."
    pip install gdown
fi

# download
gdown "https://drive.google.com/uc?export=download&id=${FILE_ID}" -O "${DEST_DIR}/${FILE_NAME}"

# check if download was successful
if [[ ! -f "${DEST_DIR}/${FILE_NAME}" ]]; then
    echo "Download failed. Please check your Google Drive link."
    exit 1
fi

# make sure the downloaded file is a ZIP archive
if ! file "${DEST_DIR}/${FILE_NAME}" | grep -q "Zip archive data"; then
    echo "Downloaded file is not a ZIP archive. Something went wrong."
    exit 1
fi

echo "File has been downloaded to ${DEST_DIR}/${FILE_NAME}."

# unzip
unzip -o "${DEST_DIR}/${FILE_NAME}" -d "$DEST_DIR"

# remove the ZIP archive
rm "${DEST_DIR}/${FILE_NAME}"

echo "Main_Dataset Folder has been downloaded and extracted to ${DEST_DIR}."