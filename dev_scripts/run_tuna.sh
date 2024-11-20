pip install tuna
python -X importtime -c "import abipy" 2> abipy_import.log
tuna abipy_import.log
