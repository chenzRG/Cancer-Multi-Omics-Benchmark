# first, make sure that there is enough disk space
# then, create a new conda environment
# finally, run this script
set -v

conda install python=3.9 -y
conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch -y
conda install numpy pandas matplotlib seaborn tqdm scikit-learn -c conda-forge -y
conda install lifelines optuna -c conda-forge -y
conda install flake8 ipdb jupyterlab tensorboard -c conda-forge -y

# pip的版本更高
pip install hydra-core