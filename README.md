# Data2016_2017
There are data I used during my last year of PhD.


## git-annex

### How to

- add krakenator rsync remote
```sh
git annex initremote krakenator_rsync type=rsync rsyncurl=krakenator.imag.fr:/home/cayek/GitAnnex/Data2016_2017 encryption=none
```

- add data
```sh
git annex add file_name
```

- copy data to krakenator_rsync
```sh
git annex copy --to krakenator_rsync
```

Data are stored on Krakenator server.



For more details see: https://git-annex.branchable.com/walkthrough/
