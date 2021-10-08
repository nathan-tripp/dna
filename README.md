# dna
entirely and totally nonfunctional at present. Waiting for things to compile and screwing off due to data loss on things I need to be working on.

Notes to self:
Helix is implemented as an inline second dimension to the array and implemented as (N+Ni)*(1+i) which apparently relates to a 45 degree rotation in the imaginary plane. This effectively comes out to be (N+N). This is done so that when the imaginary plane is copied in reverse, that some amount of the bits are not available because they would extend beyond the end of the array thus putting a finite copy count on the contents of the object. Presently unimplemented.

Using a 4 byte or similar window for the right-hand strand would allow the data to be more cleanly self-contained; but I'm tired and probably not thinking right again-- the problem was the left-right parity; the idea is that I need to not have a 1:1 mapping between the two with the link between each to overlap the boundaries of the array so that it's not possible to correctly link the left and right strands at the ends of the buffer, thus leading to corruption after N many copies of the data. At 1am, its making sense that I can use a different fragment/window size and go back to using modulus with an 8/16 byte repetition, but when I did that earlier it of course didn't work. Now I can't entirely recall why other than the left-right parity issue, but it seems like it should work... This is why I just got over the backwards references and stuck strictly to a 45 degree spiral.

Telomeres are implemented presently as a static 6 byte array with the same values at human DNA: TTAGGG.



The array is (incorrectly) described in the attached images, some modifications were made for functionality purposes but basically the lower 2 bits are metadata and intended to be used as copying space, the next two bits are the right strand, the next two bits after that are the left strand and the final 2 bits- the left-most/high-order two bits are meta-data and scratch space for the left strand. The idea is that for each copy of the data that the "fork" occurs by:


(unimplemented)
1. Finding the replication origin
2. Nicking the bits at the end of the replication origin
3. Dehelixing the bits across the entire strand ??
3. Continuously copying the left bits into the outer-most scratch space bits.
4. Inverting/flipping bits 8-4 so that 4-6 become 7-8 and vice versa.
 4a. This is the copy returned to the user.
5. Moving from the 5 prime end, or m_prime-- the process is repeated except in reverse and heuristically as a window

The replibcation origin is partly implemented, not presently helixed and stored in m_origin. A sequence of nucleotides are constructed random at object instantiation time. The idea is that the AT-ratio and GT-ratio in relation to the number of repeatitions of the N-mer sequence produces a unique series of ratios that mark the beginning of the usable space.

It likely makes a lot of sense to eventually make the sizing dynamic instead of fixed at page size. I'm mostly just playing around with an idea in my head.

Input values would be encoded along the lines of:

``    for (std::size_t idx = 0; idx < 4; idx++) {
        m_data.push_back(g_bases[(val & 0xF0) >> 6].val);
        val <<= 2;
    }``
    
yielding a 4x increase in data storage, every 2 bits incurs a byte of data.



