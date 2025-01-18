# Matchers

## Content
* [BarcodeAligner](#BarcodeAligner)
* [NEditBarcodeAligner](#NEditBarcodeAligner)

### BarcodeAligner
The Barcode aligner matches a read to a variant using a barcode sequence. <br>
What you need to use the BarcodeAligner : 
1. The design.csv file should contain a barcode column, where each varinat is associated with the corresponding variant.
1. The config.json file should contain the start and end position of the barcode.

Read will be matched to a variant only if there is a perfect match between the barcode area in the read
and the variant barcode.

 ### NEditBarcodeAligner
 The NEditBarcodeAligner matches a read to a varint using a barcode sequence.
 What you need to use the BarcodeAligner : 
1. The design.csv file should contain a barcode column, where each varinat is associated with the corresponding variant.
1. The config.json file should contain the start and end position of the barcode. as well as a "barcode_tolerance" value 
to indicate how much edit distance is allowed.

Read will be matched to a variant if there is a match between the barcode area in the read
and the variant barcode in the tolerance range, and there is only a single match.
