######## 並列計算を使えるように #########
from multiprocessing import Pool

##### 並列計算させる関数(処理):引数1つ ###
##### この場合は，引数の二乗を返す関数 ###
def nijou(x):
    print( x*x )

###### 並列計算させてみる #########
if __name__ == "__main__":
    p = Pool(4)
    p.map( nijou, range(10) )#nijouに0,1,..のそれぞれを与えて並列演算
