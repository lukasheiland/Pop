There are multiple versions of the model which were saved different times in the development cycle

	1. "ba" the standard model
	2. "ba_test" the standard model with more comments and testing options. (As of 2020–05-10 this includes a random factor m.)
	3. "ba_env" like ba, but with loc-wise parameters
	
	4. "ba-rag", a **more ragged** version compared to the original "ba"
	5. "ba-rag-ranef", like "ba-rag" but with hierarchical fitting of clusterwise parameters
	6. "ba-rect" entirely **unragged** model for simple testing purposes, otherwise entirely like ba

"ba_rag" (2) and "ba-rag-ranef" (3) are not tested in this version!
For details on  models, see file "Notes_ba.md".
