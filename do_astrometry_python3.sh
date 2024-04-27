start=`date +%s`
python full_gaia_query.py $1
python epoch_conversion_full_eq_obs.py $1
python final_transformations.py $1
python update_wcs_server.py $1
rm -rf ../extras/*.txt ../extras/*.dat ../extras/chunk*.fits ../extras/new*.fits ../extras/chunk*.pdf
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
