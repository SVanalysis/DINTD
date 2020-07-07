import dintd


if __name__ == "__main__":
    # params list
    binSize = 2000
    alpha = 0.25
    for i in range(1,2):
        try:
            reference = "chr21.fa"
            bam = "example.bam"
            str_cigar = "100M"
            binresult = "binresult"+str(binSize)+str(alpha)
            middle_result = "middle_result"+str(binSize)+str(alpha)
            
            final_result = "final_result"+str(binSize)+str(alpha)
            params = (binSize,alpha,reference,bam,str_cigar,binresult,middle_result,final_result)
            dintd.main(params)
        except OSError:
            pass
        continue
