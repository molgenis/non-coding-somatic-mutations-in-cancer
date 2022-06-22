import math
from scipy.stats.contingency import relative_risk


def calculate_relative_risk(counts_breast, counts_nonbreast, donors_breast, donors_nonbreast):
    """

    :param : 
    :param :  
    :param :        
    :return:    
    """
    RR = (counts_breast / donors_breast) / (counts_nonbreast / donors_nonbreast)
    print(RR)
    log_RR = math.log(counts_breast) - math.log(donors_breast) - math.log(counts_nonbreast) + math.log(donors_nonbreast)
    print(log_RR)
    SE_log_RR  = math.sqrt(((donors_breast-counts_breast) / (counts_breast * donors_breast)) + ((donors_nonbreast - counts_nonbreast) / (counts_nonbreast * donors_nonbreast)))
    print(SE_log_RR)
    z_a = log_RR / SE_log_RR
    print(z_a)

def main():
    calculate_relative_risk(20, 40, 286, 1952)
    print()
    result = relative_risk(20, 286, 40, 1952)
    print(result.relative_risk)
    print(result.confidence_interval(confidence_level=0.95))
    print(result.relative_risk - result.confidence_interval(confidence_level=0.95)[0])



    
if __name__ == '__main__':
    main()
