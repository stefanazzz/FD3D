subroutine inizia
include 'triffy.dec'

        v1TX=0
        v2TX=0
        v3TX=0
        s11TX=0
        s22TX=0
        s33TX=0
        s12TX=0
        s13TX=0
        s23TX=0
!=0
        v1oTX=0
        v2oTX=0
        v3oTX=0
        s11oTX=0
        s22oTX=0
        s33oTX=0
        s12oTX=0
        s13oTX=0
        s23oTX=0
        v1pTX=0
        v2pTX=0
        v3pTX=0
        s11pTX=0
        s22pTX=0
        s33pTX=0
        s12pTX=0
        s13pTX=0
        s23pTX=0
!=0
       v1TY=0
       v2TY=0
       v3TY=0
       s11TY=0
       s22TY=0
       s33TY=0
       s12TY=0
       s13TY=0
       s23TY=0
       v1oTY=0
       v2oTY=0
       v3oTY=0
       s11oTY=0
       s22oTY=0
       s33oTY=0
       s12oTY=0
       s13oTY=0
       s23oTY=0
       v1pTY=0
       v2pTY=0
       v3pTY=0
       s11pTY=0
       s22pTY=0
       s33pTY=0
       s12pTY=0
       s13pTY=0
       s23pTY=0
       v1TZ=0
       v2TZ=0
       v3TZ=0
       s11TZ=0
       s22TZ=0
       s33TZ=0
       s12TZ=0
       s13TZ=0
       s23TZ=0
       v1oTZ=0
       v2oTZ=0
       v3oTZ=0
       s11oTZ=0
       s22oTZ=0
       s33oTZ=0
       s12oTZ=0
       s13oTZ=0
       s23oTZ=0
       v1pTZ=0
       v2pTZ=0
       v3pTZ=0
       s11pTZ=0
       s22pTZ=0
       s33pTZ=0
       s12pTZ=0
       s13pTZ=0
       s23pTZ=0


        v1BX=0
        v2BX=0
        v3BX=0
        s11BX=0
        s22BX=0
        s33BX=0
        s12BX=0
        s13BX=0
        s23BX=0
!=0
        v1oBX=0
        v2oBX=0
        v3oBX=0
        s11oBX=0
        s22oBX=0
        s33oBX=0
        s12oBX=0
        s13oBX=0
        s23oBX=0
        v1pBX=0
        v2pBX=0
        v3pBX=0
        s11pBX=0
        s22pBX=0
        s33pBX=0
        s12pBX=0
        s13pBX=0
        s23pBX=0
!=0
       v1BY=0
       v2BY=0
       v3BY=0
       s11BY=0
       s22BY=0
       s33BY=0
       s12BY=0
       s13BY=0
       s23BY=0
       v1oBY=0
       v2oBY=0
       v3oBY=0
       s11oBY=0
       s22oBY=0
       s33oBY=0
       s12oBY=0
       s13oBY=0
       s23oBY=0
       v1pBY=0
       v2pBY=0
       v3pBY=0
       s11pBY=0
       s22pBY=0
       s33pBY=0
       s12pBY=0
       s13pBY=0
       s23pBY=0


       v1BZ=0
       v2BZ=0
       v3BZ=0
       s11BZ=0
       s22BZ=0
       s33BZ=0
       s12BZ=0
       s13BZ=0
       s23BZ=0
       v1oBZ=0
       v2oBZ=0
       v3oBZ=0
       s11oBZ=0
       s22oBZ=0
       s33oBZ=0
       s12oBZ=0
       s13oBZ=0
       s23oBZ=0
       v1pBZ=0
       v2pBZ=0
       v3pBZ=0
       s11pBZ=0
       s22pBZ=0
       s33pBZ=0
       s12pBZ=0
       s13pBZ=0
       s23pBZ=0
       v1TXTY=0
       v2TXTY=0
       v3TXTY=0
       s11TXTY=0
       s22TXTY=0
       s33TXTY=0
       s12TXTY=0
       s13TXTY=0
       s23TXTY=0
       v1xTXTY=0
       v2xTXTY=0
       v3xTXTY=0
       s11xTXTY=0
       s22xTXTY=0
       s33xTXTY=0
       s12xTXTY=0
       s13xTXTY=0
       v1yTXTY=0
       v2yTXTY=0
       v3yTXTY=0
       s11yTXTY=0
       s22yTXTY=0
       s33yTXTY=0
       s12yTXTY=0
       s23yTXTY=0
       v1zTXTY=0
       v2zTXTY=0
       v3zTXTY=0
       s11zTXTY=0
       s22zTXTY=0
       s33zTXTY=0
       s13zTXTY=0
       s23zTXTY=0
       v1TXBY=0
       v2TXBY=0
       v3TXBY=0
       s11TXBY=0
       s22TXBY=0
       s33TXBY=0
       s12TXBY=0
       s13TXBY=0
       s23TXBY=0
       v1xTXBY=0
       v2xTXBY=0
       v3xTXBY=0
       s11xTXBY=0
       s22xTXBY=0
       s33xTXBY=0
       s12xTXBY=0
       s13xTXBY=0
       v1yTXBY=0
       v2yTXBY=0
       v3yTXBY=0
       s11yTXBY=0
       s22yTXBY=0
       s33yTXBY=0
       s12yTXBY=0
       s23yTXBY=0
       v1zTXBY=0
       v2zTXBY=0
       v3zTXBY=0
       s11zTXBY=0
       s22zTXBY=0
       s33zTXBY=0
       s13zTXBY=0
       s23zTXBY=0

       v1BXTY=0
       v2BXTY=0
       v3BXTY=0
       s11BXTY=0
       s22BXTY=0
       s33BXTY=0
       s12BXTY=0
       s13BXTY=0
       s23BXTY=0
       v1xBXTY=0
       v2xBXTY=0
       v3xBXTY=0
       s11xBXTY=0
       s22xBXTY=0
       s33xBXTY=0
       s12xBXTY=0
       s13xBXTY=0
       v1yBXTY=0
       v2yBXTY=0
       v3yBXTY=0
       s11yBXTY=0
       s22yBXTY=0
       s33yBXTY=0
       s12yBXTY=0
       s23yBXTY=0
       v1zBXTY=0
       v2zBXTY=0
       v3zBXTY=0
       s11zBXTY=0
       s22zBXTY=0
       s33zBXTY=0
       s13zBXTY=0
       s23zBXTY=0

       v1BXBY=0
       v2BXBY=0
      v3BXBY=0
       s11BXBY=0
       s22BXBY=0
       s33BXBY=0
       s12BXBY=0
       s13BXBY=0
       s23BXBY=0
       v1xBXBY=0
       v2xBXBY=0
       v3xBXBY=0
       s11xBXBY=0
       s22xBXBY=0
       s33xBXBY=0
       s12xBXBY=0
       s13xBXBY=0
       v1yBXBY=0
       v2yBXBY=0
       v3yBXBY=0
       s11yBXBY=0
       s22yBXBY=0
       s33yBXBY=0
       s12yBXBY=0
       s23yBXBY=0
       v1zBXBY=0
       v2zBXBY=0
       v3zBXBY=0
       s11zBXBY=0
       s22zBXBY=0
       s33zBXBY=0
       s13zBXBY=0
       s23zBXBY=0

       v1TXTZ=0
       v2TXTZ=0
       v3TXTZ=0
       s11TXTZ=0
       s22TXTZ=0
       s33TXTZ=0
       s12TXTZ=0
       s13TXTZ=0
       v1TXTZ=0
       v2TXTZ=0
       v3TXTZ=0
       s11TXTZ=0
       s22TXTZ=0
       s33TXTZ=0
       s12TXTZ=0
       s13TXTZ=0
       s23TXTZ=0
       v1xTXTZ=0
       v2xTXTZ=0
       v3xTXTZ=0
       s11xTXTZ=0
       s22xTXTZ=0
       s33xTXTZ=0
       s12xTXTZ=0
       s13xTXTZ=0
       v1yTXTZ=0
       v2yTXTZ=0
       v3yTXTZ=0
       s11yTXTZ=0
       s22yTXTZ=0
       s33yTXTZ=0
       s12yTXTZ=0
       s23yTXTZ=0
       v1zTXTZ=0
       v2zTXTZ=0
       v3zTXTZ=0
       s11zTXTZ=0
       s22zTXTZ=0
       s33zTXTZ=0
       s13zTXTZ=0
       s23zTXTZ=0

       v1TXBZ=0
       v2TXBZ=0
       v3TXBZ=0
       s11TXBZ=0
      s22TXBZ=0
       s33TXBZ=0
       s12TXBZ=0
       s13TXBZ=0
       v1TXBZ=0
       v2TXBZ=0
       v3TXBZ=0
       s11TXBZ=0
       s22TXBZ=0
       s33TXBZ=0
       s12TXBZ=0
       s13TXBZ=0
       s23TXBZ=0
       v1xTXBZ=0
       v2xTXBZ=0
       v3xTXBZ=0
       s11xTXBZ=0
       s22xTXBZ=0
       s33xTXBZ=0
       s12xTXBZ=0
       s13xTXBZ=0
       v1yTXBZ=0
       v2yTXBZ=0
       v3yTXBZ=0
       s11yTXBZ=0
       s22yTXBZ=0
       s33yTXBZ=0
       s12yTXBZ=0
       s23yTXBZ=0
       v1zTXBZ=0
       v2zTXBZ=0
       v3zTXBZ=0
       s11zTXBZ=0
       s22zTXBZ=0
       s33zTXBZ=0
       s13zTXBZ=0
       s23zTXBZ=0
       v1BXTZ=0
       v2BXTZ=0
       v3BXTZ=0
       s11BXTZ=0
       s22BXTZ=0
       s33BXTZ=0
       s12BXTZ=0
       s13BXTZ=0
       v1BXTZ=0
       v2BXTZ=0
       v3BXTZ=0
       s11BXTZ=0
       s22BXTZ=0
       s33BXTZ=0
       s12BXTZ=0
       s13BXTZ=0
       s23BXTZ=0
       v1xBXTZ=0
       v2xBXTZ=0
       v3xBXTZ=0
       s11xBXTZ=0
       s22xBXTZ=0
       s33xBXTZ=0
       s12xBXTZ=0
       s13xBXTZ=0
       v1yBXTZ=0
       v2yBXTZ=0
       v3yBXTZ=0
       s11yBXTZ=0
       s22yBXTZ=0
       s33yBXTZ=0
       s12yBXTZ=0
       s23yBXTZ=0
       v1zBXTZ=0
       v2zBXTZ=0
       v3zBXTZ=0
       s11zBXTZ=0
       s22zBXTZ=0
       s33zBXTZ=0
       s13zBXTZ=0
       s23zBXTZ=0
       v1BXBZ=0
       v2BXBZ=0
      v3BXBZ=0
       s11BXBZ=0
       s22BXBZ=0
       s33BXBZ=0
       s12BXBZ=0
       s13BXBZ=0
       v1BXBZ=0
       v2BXBZ=0
       v3BXBZ=0
       s11BXBZ=0
       s22BXBZ=0
       s33BXBZ=0
       s12BXBZ=0
       s13BXBZ=0
       s23BXBZ=0
       v1xBXBZ=0
       v2xBXBZ=0
       v3xBXBZ=0
       s11xBXBZ=0
       s22xBXBZ=0
       s33xBXBZ=0
       s12xBXBZ=0
       s13xBXBZ=0
       v1yBXBZ=0
       v2yBXBZ=0
       v3yBXBZ=0
       s11yBXBZ=0
       s22yBXBZ=0
       s33yBXBZ=0
       s12yBXBZ=0
       s23yBXBZ=0
       v1zBXBZ=0
       v2zBXBZ=0
       v3zBXBZ=0
       s11zBXBZ=0
       s22zBXBZ=0
       s33zBXBZ=0
       s13zBXBZ=0
       s23zBXBZ=0

       v1TYTZ=0
       v2TYTZ=0
       v3TYTZ=0
       s11TYTZ=0
       s22TYTZ=0
       s33TYTZ=0
       s12TYTZ=0
       s13TYTZ=0
       s23TYTZ=0
       v1xTYTZ=0
       v2xTYTZ=0
       v3xTYTZ=0
       s11xTYTZ=0
       s22xTYTZ=0
       s33xTYTZ=0
       s12xTYTZ=0
       s13xTYTZ=0
       v1yTYTZ=0
       v2yTYTZ=0
       v3yTYTZ=0
       s11yTYTZ=0
       s22yTYTZ=0
       s33yTYTZ=0
       s12yTYTZ=0
       s23yTYTZ=0
       v1zTYTZ=0
       v2zTYTZ=0
       v3zTYTZ=0
       s11zTYTZ=0
       s22zTYTZ=0
       s33zTYTZ=0
       s13zTYTZ=0
       s23zTYTZ=0

       v1TYBZ=0
       v2TYBZ=0
       v3TYBZ=0
       s11TYBZ=0
       s22TYBZ=0
       s33TYBZ=0
       s12TYBZ=0
       s13TYBZ=0
       s23TYBZ=0
       v1xTYBZ=0
       v2xTYBZ=0
       v3xTYBZ=0
       s11xTYBZ=0
       s22xTYBZ=0
       s33xTYBZ=0
       s12xTYBZ=0
       s13xTYBZ=0
       v1yTYBZ=0
       v2yTYBZ=0
       v3yTYBZ=0
       s11yTYBZ=0
       s22yTYBZ=0
       s33yTYBZ=0
       s12yTYBZ=0
       s23yTYBZ=0
       v1zTYBZ=0
       v2zTYBZ=0
       v3zTYBZ=0
       s11zTYBZ=0
       s22zTYBZ=0
       s33zTYBZ=0
       s13zTYBZ=0
       s23zTYBZ=0

       v1BYTZ=0
       v2BYTZ=0
       v3BYTZ=0
       s11BYTZ=0
      s22BYTZ=0
       s33BYTZ=0
       s12BYTZ=0
       s13BYTZ=0
       s23BYTZ=0
       v1xBYTZ=0
       v2xBYTZ=0
       v3xBYTZ=0
       s11xBYTZ=0
       s22xBYTZ=0
       s33xBYTZ=0
       s12xBYTZ=0
       s13xBYTZ=0
       v1yBYTZ=0
       v2yBYTZ=0
       v3yBYTZ=0
       s11yBYTZ=0
       s22yBYTZ=0
       s33yBYTZ=0
       s12yBYTZ=0
       s23yBYTZ=0
       v1zBYTZ=0
       v2zBYTZ=0
       v3zBYTZ=0
       s11zBYTZ=0
       s22zBYTZ=0
       s33zBYTZ=0
       s13zBYTZ=0
       s23zBYTZ=0

       v1BYBZ=0
       v2BYBZ=0
       v3BYBZ=0
       s11BYBZ=0
       s22BYBZ=0
       s33BYBZ=0
       s12BYBZ=0
       s13BYBZ=0
       s23BYBZ=0
       v1xBYBZ=0
       v2xBYBZ=0
       v3xBYBZ=0
       s11xBYBZ=0
       s22xBYBZ=0
       s33xBYBZ=0
       s12xBYBZ=0
       s13xBYBZ=0
       v1yBYBZ=0
       v2yBYBZ=0
       v3yBYBZ=0
       s11yBYBZ=0
       s22yBYBZ=0
       s33yBYBZ=0
       s12yBYBZ=0
       s23yBYBZ=0
       v1zBYBZ=0
       v2zBYBZ=0
       v3zBYBZ=0
       s11zBYBZ=0
       s22zBYBZ=0
       s33zBYBZ=0
       s13zBYBZ=0
       s23zBYBZ=0

       v1TXTYTZ=0
       v2TXTYTZ=0
       v3TXTYTZ=0
       s11TXTYTZ=0
       s22TXTYTZ=0
       s33TXTYTZ=0
       s12TXTYTZ=0
       s13TXTYTZ=0
       s23TXTYTZ=0
       v1xTXTYTZ=0
       v2xTXTYTZ=0
       v3xTXTYTZ=0
       s11xTXTYTZ=0
       s22xTXTYTZ=0
       s33xTXTYTZ=0
       s12xTXTYTZ=0
       s13xTXTYTZ=0
       v1yTXTYTZ=0
       v2yTXTYTZ=0
       v3yTXTYTZ=0
       s11yTXTYTZ=0
       s22yTXTYTZ=0
       s33yTXTYTZ=0
       s12yTXTYTZ=0
       s23yTXTYTZ=0
       v1zTXTYTZ=0
       v2zTXTYTZ=0
       v3zTXTYTZ=0
       s11zTXTYTZ=0
       s22zTXTYTZ=0
       s33zTXTYTZ=0
       s13zTXTYTZ=0
       s23zTXTYTZ=0

       v1TXTYBZ=0
       v2TXTYBZ=0
       v3TXTYBZ=0
       s11TXTYBZ=0
       s22TXTYBZ=0
       s33TXTYBZ=0
       s12TXTYBZ=0
       s13TXTYBZ=0
       s23TXTYBZ=0
       v1xTXTYBZ=0
       v2xTXTYBZ=0
       v3xTXTYBZ=0
       s11xTXTYBZ=0
       s22xTXTYBZ=0
       s33xTXTYBZ=0
       s12xTXTYBZ=0
       s13xTXTYBZ=0
       v1yTXTYBZ=0
       v2yTXTYBZ=0
       v3yTXTYBZ=0
       s11yTXTYBZ=0
       s22yTXTYBZ=0
       s33yTXTYBZ=0
       s12yTXTYBZ=0
       s23yTXTYBZ=0
       v1zTXTYBZ=0
       v2zTXTYBZ=0
       v3zTXTYBZ=0
       s11zTXTYBZ=0
       s22zTXTYBZ=0
       s33zTXTYBZ=0
       s13zTXTYBZ=0
       s23zTXTYBZ=0
       v1TXBYTZ=0
       v2TXBYTZ=0
       v3TXBYTZ=0
       s11TXBYTZ=0
       s22TXBYTZ=0
       s33TXBYTZ=0
       s12TXBYTZ=0
       s13TXBYTZ=0
       s23TXBYTZ=0
       v1xTXBYTZ=0
       v2xTXBYTZ=0
       v3xTXBYTZ=0
       s11xTXBYTZ=0
       s22xTXBYTZ=0
       s33xTXBYTZ=0
       s12xTXBYTZ=0
       s13xTXBYTZ=0
       v1yTXBYTZ=0
       v2yTXBYTZ=0
       v3yTXBYTZ=0
       s11yTXBYTZ=0
       s22yTXBYTZ=0
       s33yTXBYTZ=0
       s12yTXBYTZ=0
       s23yTXBYTZ=0
       v1zTXBYTZ=0
       v2zTXBYTZ=0
       v3zTXBYTZ=0
       s11zTXBYTZ=0
       s22zTXBYTZ=0
       s33zTXBYTZ=0
       s13zTXBYTZ=0
       s23zTXBYTZ=0
       v1TXBYBZ=0
       v2TXBYBZ=0
       v3TXBYBZ=0
       s11TXBYBZ=0
       s22TXBYBZ=0
       s33TXBYBZ=0
       s12TXBYBZ=0
       s13TXBYBZ=0
       s23TXBYBZ=0
       v1xTXBYBZ=0
       v2xTXBYBZ=0
       v3xTXBYBZ=0
       s11xTXBYBZ=0
       s22xTXBYBZ=0
       s33xTXBYBZ=0
       s12xTXBYBZ=0
       s13xTXBYBZ=0
       v1yTXBYBZ=0
       v2yTXBYBZ=0
       v3yTXBYBZ=0
       s11yTXBYBZ=0
       s22yTXBYBZ=0
       s33yTXBYBZ=0
       s12yTXBYBZ=0
       s23yTXBYBZ=0
       v1zTXBYBZ=0
       v2zTXBYBZ=0
       v3zTXBYBZ=0
       s11zTXBYBZ=0
       s22zTXBYBZ=0
       s33zTXBYBZ=0
       s13zTXBYBZ=0
       s23zTXBYBZ=0
       v1BXTYTZ=0
       v2BXTYTZ=0
       v3BXTYTZ=0
       s11BXTYTZ=0
       s22BXTYTZ=0
       s33BXTYTZ=0
       s12BXTYTZ=0
       s13BXTYTZ=0
       s23BXTYTZ=0
       v1xBXTYTZ=0
       v2xBXTYTZ=0
       v3xBXTYTZ=0
       s11xBXTYTZ=0
       s22xBXTYTZ=0
       s33xBXTYTZ=0
       s12xBXTYTZ=0
       s13xBXTYTZ=0
       v1yBXTYTZ=0
       v2yBXTYTZ=0
       v3yBXTYTZ=0
       s11yBXTYTZ=0
       s22yBXTYTZ=0
       s33yBXTYTZ=0
       s12yBXTYTZ=0
       s23yBXTYTZ=0
       v1zBXTYTZ=0
       v2zBXTYTZ=0
       v3zBXTYTZ=0
       s11zBXTYTZ=0
       s22zBXTYTZ=0
       s33zBXTYTZ=0
       s13zBXTYTZ=0
       s23zBXTYTZ=0
       v1BXTYBZ=0
       v2BXTYBZ=0
       v3BXTYBZ=0
       s11BXTYBZ=0
       s22BXTYBZ=0
       s33BXTYBZ=0
       s12BXTYBZ=0
       s13BXTYBZ=0
       s23BXTYBZ=0
       v1xBXTYBZ=0
       v2xBXTYBZ=0
       v3xBXTYBZ=0
       s11xBXTYBZ=0
       s22xBXTYBZ=0
       s33xBXTYBZ=0
       s12xBXTYBZ=0
       s13xBXTYBZ=0
       v1yBXTYBZ=0
       v2yBXTYBZ=0
       v3yBXTYBZ=0
       s11yBXTYBZ=0
       s22yBXTYBZ=0
       s33yBXTYBZ=0
       s12yBXTYBZ=0
       s23yBXTYBZ=0
       v1zBXTYBZ=0
       v2zBXTYBZ=0
       v3zBXTYBZ=0
       s11zBXTYBZ=0
       s22zBXTYBZ=0
       s33zBXTYBZ=0
       s13zBXTYBZ=0
       s23zBXTYBZ=0
       v1BXBYTZ=0
       v2BXBYTZ=0
       v3BXBYTZ=0
       s11BXBYTZ=0
       s22BXBYTZ=0
       s33BXBYTZ=0
       s12BXBYTZ=0
       s13BXBYTZ=0
       s23BXBYTZ=0
       v1xBXBYTZ=0
       v2xBXBYTZ=0
       v3xBXBYTZ=0
       s11xBXBYTZ=0
       s22xBXBYTZ=0
       s33xBXBYTZ=0
       s12xBXBYTZ=0
       s13xBXBYTZ=0
       v1yBXBYTZ=0
       v2yBXBYTZ=0
       v3yBXBYTZ=0
       s11yBXBYTZ=0
       s22yBXBYTZ=0
       s33yBXBYTZ=0
       s12yBXBYTZ=0
       s23yBXBYTZ=0
       v1zBXBYTZ=0
       v2zBXBYTZ=0
      v3zBXBYTZ=0
       s11zBXBYTZ=0
       s22zBXBYTZ=0
       s33zBXBYTZ=0
       s13zBXBYTZ=0
       s23zBXBYTZ=0
       v1BXBYBZ=0
       v2BXBYBZ=0
       v3BXBYBZ=0
       s11BXBYBZ=0
       s22BXBYBZ=0
       s33BXBYBZ=0
       s12BXBYBZ=0
       s13BXBYBZ=0
       s23BXBYBZ=0
       v1xBXBYBZ=0
       v2xBXBYBZ=0
       v3xBXBYBZ=0
       s11xBXBYBZ=0
       s22xBXBYBZ=0
       s33xBXBYBZ=0
       s12xBXBYBZ=0
       s13xBXBYBZ=0
       v1yBXBYBZ=0
       v2yBXBYBZ=0
       v3yBXBYBZ=0
       s11yBXBYBZ=0
       s22yBXBYBZ=0
       s33yBXBYBZ=0
       s12yBXBYBZ=0
       s23yBXBYBZ=0
       v1zBXBYBZ=0
       v2zBXBYBZ=0
       v3zBXBYBZ=0
       s11zBXBYBZ=0
       s22zBXBYBZ=0
       s33zBXBYBZ=0
       s13zBXBYBZ=0
       s23zBXBYBZ=0
return
end
