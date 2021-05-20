### Upload Files to AWS S3

From crosman.genetics.utah.edu Dec 18 2018.

```
ssh bmoore@crosman.genetics.utah.edu
cd /archive01/bmoore/Projects/18-12-01_Build_Gnomad_VVP_Bckgnd
ls Gnomad_2.1_VVP_Bckgnd_controls_vvp_bkgd.* | nohup parallel 'aws s3 cp --profile omicia {} s3://vvp-gnomad-r2.1-background 2> aws-s3-cp_{}.error' &
```
# catherpes
