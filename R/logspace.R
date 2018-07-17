logspace<-function (a, b, n = 50, base = 10) 
{
  if (b == pi && base == 10) {
    b <- log(b, base = base)
  }
  base^seq(from = a, to = b, length.out = n)
}