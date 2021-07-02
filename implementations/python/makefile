test:
	pytest -r a -v tests --cov --cov-report=html --cov-report term

rebuild_test_bundle:
	mzspeclib convert -f text tests/test_data/chinese_hamster_hcd_selected_head.msp tests/test_data/chinese_hamster_hcd_selected_head.mzlb.txt
	mzspeclib convert -f json tests/test_data/chinese_hamster_hcd_selected_head.msp tests/test_data/chinese_hamster_hcd_selected_head.mzlb.json
	mzspeclib convert -f json "tests/test_data/complex_interpretations_with_members.mzlb.txt" tests/test_data/complex_interpretations_with_members.mzlb.json
	python tests/test_data/generate_annotations.py