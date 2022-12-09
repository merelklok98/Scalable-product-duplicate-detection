Product Duplicate Detection based on Product Titles and Features using LSH and Jaccard similarity

An algorithm for finding duplicate products in a dataset of TVs of four different Webshops. Users can more easily compare various items by using this algorithm. Applying LSH to identify potential candidate pairings reduces computation time. Jaccard Similarity is then used to assess the accuracy of these pairs and define the final predicted duplicates.

The data is cleaned up in the first section of the code. The features that will be utilized for the LSH stage are then chosen out of the titles. In this instance, only the words with alphabetic and numeric characters are chosen. With min-hashing the dimension of the created binary matrix is reduced and with the LSH method a candidate pairs-matrix is created. In the final step the Jaccard dissimilarity of the candidate pairs is determined to point out the real predicted duplicates. 

The features chosen, the cleaning process, the calling or not of the built functions, and the threshold values can all be changed by the user. The final results are evaluated with the F1-score to determine how well the detection method performed.

This code has been written on behalf of Erasmus University by Merel Klok 
