# Data2016_2017
There are data I used during my last year of PhD.


## git-annex

In new git repo (git bare repo too)
```
git annex init
```
Data are on TIMC-BCM-15 and krakenator.


### How to

- add krakenator remote
```sh
git remote add krakenator krakenator.imag.fr:/home/cayek/Project/Data2016_2017
```

- add data
```
git annex add file_name
```

- copy data to krakenator
```
git annex copy --to krakenator
```

- get data from krakenator_rsync
```
git annex sync krakenator
git annex get file_name
```



For more details see: https://git-annex.branchable.com/walkthrough/
