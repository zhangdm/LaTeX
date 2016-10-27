# # =========== grammer ====# #
dunif(x, min = 0, max = 1, log = FALSE)

punif(q, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)  

qunif(p, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)

runif(n, min = 0, max = 1) 

# # =========example ==========# #
u <- runif(20)

## The following relations always hold :
punif(u) == u
dunif(u) == 1
var(runif(10000))  #- ~ = 1/12 = .08333