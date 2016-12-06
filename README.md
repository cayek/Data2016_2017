# Data2016_2017
There are data I used during my last year of PhD.


## git-annex

```
git annex init
```

### How to

- add krakenator remote
```sh
git remote add krakenator krakenator.imag.fr:/home/cayek/GitRepo/Data2016_2017.git
```

- add data
```
git annex add file_name
```

- copy data to krakenator_rsync
```
git annex copy --to krakenator_rsync
```

- get data from krakenator_rsync
```
git annex sync krakenator_rsync
git annex get file_name
```
Data are stored on Krakenator server.



For more details see: https://git-annex.branchable.com/walkthrough/
