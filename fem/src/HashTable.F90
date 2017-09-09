!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 09 Aug 1998   (Original C-version) Genesis
! *  Original Date: 25 Oct 2000   Fortran 90/95 version
! *
! *****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!-------------------------------------------------------
!>  Hash table build & lookup routines.
!-------------------------------------------------------

MODULE HashTable

  USE Lists

  IMPLICIT NONE

  TYPE HashValue_t
     CHARACTER(LEN=MAX_NAME_LEN) :: Block,Type
  END TYPE HashValue_t


  TYPE HashEntry_t
     TYPE(HashEntry_t), POINTER  :: Next
     CHARACTER(LEN=MAX_NAME_LEN) :: Key
     TYPE(HashValue_t), POINTER  :: Value
  END TYPE HashEntry_t


  TYPE HashBucket_t
     TYPE(HashEntry_t), POINTER :: Head
  END TYPE HashBucket_t


  TYPE HashTable_t
     INTEGER :: BucketSize, CurrentBucket, TotalEntries
     INTEGER :: MaxAvgEntries
     TYPE(HashBucket_t), POINTER :: Bucket(:)
     TYPE(HashEntry_t),  POINTER :: CurrentEntry
  END TYPE HashTable_t


CONTAINS

!-----------------------------------------------------------------------
! 
!  Call: TYPE(HashTable_t), POINTER hash = HashCreate( &
!       InitialBucketSize, MaxAvgEntries, EqualKeys )
! 
!> Initialize a hash table given initial bucket size. The size of the
!> bucket is rounded up to next power of two. The bucket doubles in size
!> whenever the size of the hash table grows over "MaxAvgEntries"
!> entries / bucket on the average. Keep the "MaxAvgEntries" small enough
!> (ordinarily from 3 entries up ?) to keep the hash table build & lookup
!> reasonably quick.
! 
!> The hash table entries contain a "key" and an associated "value".
!> Currently the value structures are not copied to the hash table,
!> but only pointers to them are being stored. So one should not alter the
!> memory where the components are actually stored after adding them to the
!> hash table. Note that, potential problems in this respect are,
!> among others, automatic variables in functions. (C-Version: Bright side
!> is, that the "value" entries may be whatsoever as you manage them
!> yourself...)
! 
!  (Again, the C-version:
!  The third and fourth arguments may be used to give user routines to
!  compare equality of two keys and provide an index to the bucket array
!  given key and array size respectively (these are the only place where
!  any reference to the format of the keys is made). If given as NULL
!  pointers, the default string comparison and hash functions routines
!  respectively are implied. )
! 
!-----------------------------------------------------------------------
  FUNCTION HashCreate( InitialBucketSize, MaxAvgEntries ) RESULT(hash)
    TYPE(HashTable_t), POINTER :: Hash
    INTEGER :: InitialBucketSize, MaxAvgEntries

    INTEGER :: i, RoundBits, Stat

    NULLIFY( Hash )
    IF ( InitialBucketSize <= 0 ) THEN
       WRITE( Message, * )  'HashCreate: invalid initial size given: ',InitialBucketSize
       CALL Error( 'HashCreate', Message )
       RETURN
    END IF

!   /*
!    *  Round bucket size up to next largest power of two...
!    */
    RoundBits = CEILING( LOG(1.0d0*InitialBucketSize) / LOG(2.0d0) )
 
!   /*
!    * Allocate the table and initialize the table entries...
!    */
    ALLOCATE( Hash )

    Hash % BucketSize = 2**RoundBits;

!   /*
!    *  Allocate the bucket array
!    */
    ALLOCATE( Hash % Bucket( Hash % BucketSize ), STAT=stat )

    IF ( stat /= 0 ) THEN
       CALL Error( 'HashCreate', &
               'Hash table initialize error: unable to allocate bucket.' )
       DEALLOCATE( Hash )
       NULLIFY( Hash )

       RETURN
    END IF

    DO i=1,Hash % BucketSize
       NULLIFY( Hash % Bucket(i) % Head )
    END DO

    Hash % TotalEntries  = 0
    Hash % MaxAvgEntries = MaxAvgEntries;

!   /*
!   *  the key comparison routine
!    */
!   if ( EqualKeys )
!     hash->EqualKeys = EqualKeys;
!   else
!     hash->EqualKeys = (int (*)(void *,void *))HashEqualKeys;
!
!   /*
!    *  the hash table index generation routine
!    */
!   if ( HashFunc )
!     hash->HashFunc = HashFunc;
!   else
!     hash->HashFunc = (int (*)(void *,int))HashStringFunc;
   END FUNCTION HashCreate

!--------------------------------------------------------------------------
!  Call:index = HashStringFunc( key, mask )
! 
!> Generate index to a hash table from given string. Hash table size
!> is assumed to be a power of two.
!--------------------------------------------------------------------------
   FUNCTION HashStringFunc( key, size ) RESULT(Ind)
     CHARACTER(LEN=*) :: key
     INTEGER :: Size, Ind

     INTEGER :: i,keylen

     DO keylen=LEN(key),1,-1
       IF ( key(keylen:keylen) /= ' ' ) EXIT
     END DO

     Ind = 0
     DO i=1,keylen
        Ind = Ind*8 + ICHAR(key(i:i))
     END DO

     Ind = IAND( Ind, size-1 ) + 1
   END FUNCTION HashStringFunc

!--------------------------------------------------------------------------
! 
! Call: equal = HashEqualKeys( key1, key2 )
! 
!> Return equality of given two strings.
!> This is for internal use only.
!--------------------------------------------------------------------------
   FUNCTION HashEqualKeys( key1,key2 ) RESULT(equal)

     CHARACTER(LEN=*) :: key1,key2
     LOGICAL :: equal
     INTEGER :: n1,n2

     equal = .FALSE.

     n1 = LEN_TRIM(key1)
     n2 = LEN_TRIM(key2)
     IF ( n1 /= n2 ) RETURN

     equal = key1(1:n1) == key2(1:n1)
   END FUNCTION HashEqualKeys

!--------------------------------------------------------------------------
!  Call: entry = HashFind( hash, key, bucket )
! 
!> Search for a key from a hash table, return value is pointer to
!> the entry or NULL if not found. Bucket number of the entry
!> (if found) is given in int *bucket.
!> This is for internal use only.
!
!--------------------------------------------------------------------------
   FUNCTION HashFind( hash, key, n ) RESULT( entry )

     TYPE(HashTable_t), POINTER :: Hash
     TYPE(HashEntry_t), POINTER :: Entry, ptr
     INTEGER :: n
     CHARACTER(LEN=*) :: key

!    /*
!     * Get bucket number for this key
!     */
     n = HashStringFunc( key, Hash % BucketSize )

!    /*
!     * Now go through bucket entries, and check if it is to be found
!     */
    ptr => Hash % Bucket(n) % Head
    DO WHILE( ASSOCIATED(ptr) )
       IF ( HashEqualKeys( key, ptr % key ) ) EXIT
       ptr => ptr % Next
    END DO
    Entry => ptr
  END FUNCTION HashFind

!--------------------------------------------------------------------------
! Call:  HashAdd( HashTable_t *hash, void *key,void *value )
! 
!> Add an entry to a hash table. If the key is already in the table
!> just change the "value" pointer.
! 
!> The hash table entries contain a "key" and an associated "value".
!> Currently the (key) and value entries are not copied to the hash table,
!> but only pointers to them are being stored. So one should not alter the
!> memory where the components are actually stored after adding them to the
!> hash table. Note that, potential problems in this respect are, for
!> example, automatic variables in functions.
!>
!  Return value is success or not...
! 
!--------------------------------------------------------------------------
  RECURSIVE FUNCTION HashAdd( hash, key, value ) RESULT(Success)
    TYPE(HashTable_t), POINTER :: hash
    CHARACTER(LEN=*) :: key
    LOGICAL :: Success
    TYPE(HashValue_t), POINTER :: Value

    INTEGER :: stat
    INTEGER :: n, keylen, count = 0
    TYPE(HashEntry_t), POINTER :: Entry

    Success = .TRUE.

    entry => HashFind( hash,key, n )

    IF ( ASSOCIATED( entry ) ) THEN
!      /*
!       * Already in, change the value pointer
!       */
       Entry % Value => Value
    ELSE
!       /*
!        * not found add new...
!        */
       ALLOCATE( Entry, STAT=stat )

       IF ( stat /= 0 ) THEN
          CALL Error( 'HashAdd', 'add failed: unable to allocate ' // &
               '(a few bytes of) memory ?' )
          RETURN
       END IF

       Entry % Next  => Hash % Bucket(n) % Head
       Entry % Value => Value
       Entry % Key = ' '
       DO keylen=LEN(key),1,-1
          IF ( key(keylen:keylen) /= ' ' ) EXIT
       END DO
       Entry % Key(1:keylen) = key(1:keylen)

       Hash % Bucket(n) % Head => Entry
       Hash % TotalEntries = Hash % TotalEntries + 1
       
       IF ( Hash % TotalEntries > Hash % MaxAvgEntries*Hash % BucketSize ) THEN
          Success = HashRebuild( Hash )
       END IF
    END IF
  END FUNCTION HashAdd

!--------------------------------------------------------------------------
!  Call: HashRemove( HashTable_t *hash, void *key )
! 
!> Remove an entry from a hash table given key of the entry.
!--------------------------------------------------------------------------
  SUBROUTINE HashRemove( Hash, key )
    TYPE(HashTable_t), POINTER :: Hash
    CHARACTER(LEN=*) :: Key

    TYPE(HashEntry_t), POINTER :: entry,prev
    INTEGER :: k,n

    IF ( .NOT. ASSOCIATED(hash) ) RETURN

!    /*
!     *  get bucket number for this key
!     */
    n = HashStringFunc( key, hash % BucketSize );

!    /*
!     * Now go through bucket entries, and check if it's there
!     */
    NULLIFY( Prev )
    Entry => Hash % Bucket(n) % Head
    DO WHILE( ASSOCIATED( Entry ) )
!       /*
!        * if key in, remove
!        */
       IF ( HashEqualKeys( key,entry % key ) ) THEN
          IF ( ASSOCIATED(prev) ) THEN
            prev % next => entry % next
          ELSE
            hash % Bucket(n) % Head => Entry % next
         END IF

         DEALLOCATE(Entry)

         hash  % TotalEntries = Hash % TotalEntries - 1
         RETURN
      END IF
      Prev  => Entry
      Entry => Entry % Next
   END DO
 END SUBROUTINE HashRemove

!--------------------------------------------------------------------------
! Call: HashClean( HashTable_t *hash )
! 
!> Clean all entries from the hash table, the bucket array is kept.
!> One may start refilling the hash table directly after cleaning.
!--------------------------------------------------------------------------
 SUBROUTINE HashClean( hash )
   TYPE(HashTable_t), POINTER :: Hash

   TYPE(HashEntry_t), POINTER :: ptr,ptr1
   INTEGER :: i

   IF ( .NOT.ASSOCIATED(hash) ) RETURN

   DO i=1,hash % BucketSize
      ptr => hash % Bucket(i) % Head
      DO WHILE( ASSOCIATED(ptr) )
         ptr1 => ptr % next
         DEALLOCATE( ptr )
         ptr => ptr1
      END DO
      NULLIFY( hash % Bucket(i) % Head )
   END DO
   hash % TotalEntries  =  0
   hash % CurrentBucket =  0
   NULLIFY( hash % CurrentEntry )
 END SUBROUTINE HashClean

!--------------------------------------------------------------------------
!  Call: HashDelete( HashTable_t *hash )
!
!> Delete a hash table by removing all the entries and freeing the
!> bucket and hash structures.
!--------------------------------------------------------------------------
 SUBROUTINE HashDelete( Hash )
   TYPE(HashTable_t), POINTER :: Hash

   IF ( ASSOCIATED( Hash ) ) THEN
      IF ( ASSOCIATED( Hash % Bucket ) ) THEN
         CALL HashClean( Hash )
         DEALLOCATE( Hash % Bucket )
      END IF
      DEALLOCATE( Hash )
   END IF
 END SUBROUTINE HashDelete

!--------------------------------------------------------------------------
! Call: HashRebuild( HashTable_t *hash )
! 
!> Rebuild a hash table using a larger bucket array.
!> This is for internal use only.
!--------------------------------------------------------------------------
 RECURSIVE FUNCTION HashRebuild( hash ) RESULT(Success)
   TYPE(HashTable_t), POINTER :: Hash, NewHash
   LOGICAL :: Success

   TYPE(HashEntry_t), POINTER :: entry
   INTEGER :: i

   Success = .FALSE.
   IF ( .NOT.ASSOCIATED(Hash) ) RETURN

   NewHash => HashCreate( 4*Hash % BucketSize, Hash % MaxAvgEntries )
   IF ( .NOT. ASSOCIATED( Newhash ) ) RETURN

   DO i=1,Hash % BucketSize
      Entry => Hash % Bucket(i) % Head
      DO WHILE( ASSOCIATED( Entry ) )
         IF ( .NOT. HashAdd( Newhash, Entry % Key, Entry % Value ) ) RETURN
         Entry => Entry % Next
      END DO
   END DO

   CALL HashClean( Hash )
   DEALLOCATE( Hash % Bucket )

   DEALLOCATE( Hash )
   Hash => NewHash
   Success = .TRUE.
 END FUNCTION HashRebuild

!--------------------------------------------------------------------------
! Call: void *value = HashValue( HashTable_t *hash, void *key )
! 
!> Given a "key" to hash table return pointer to the "value" memory or
!> NULL if not found in the table.
!--------------------------------------------------------------------------
 FUNCTION HashValue( Hash, key ) RESULT(Value)
   TYPE(HashTable_t), POINTER :: Hash
   CHARACTER(LEN=*) :: Key
   TYPE(HashValue_t), POINTER :: Value

   INTEGER ::  n
   TYPE(HashEntry_t), POINTER :: Entry

   NULLIFY( Value )
   Entry => HashFind( Hash, key, n )
   IF ( ASSOCIATED( Entry ) ) Value => Entry % Value
 END FUNCTION HashValue

!--------------------------------------------------------------------------
! Call: HashInitWalk( HashTable_t *hash )
!
!> Initialize hash table walk through.
!--------------------------------------------------------------------------
 SUBROUTINE HashInitWalk( Hash )
   TYPE(HashTable_t), POINTER :: Hash
   Hash % CurrentBucket = 0
   NULLIFY( Hash % CurrentEntry )
 END SUBROUTINE HashInitWalk

!--------------------------------------------------------------------------
!! Call: HashEntry_t *entry = HashNext( HashTable_t *hash )
!!
!> Return pointer to "next" entry in a hash table. The walk must be
!> initialized with a call to HashInitWalk. The "key" and "value" of
!> the table entry may be referenced as follows:
!--------------------------------------------------------------------------
 FUNCTION HashNext( Hash ) RESULT(Entry)
   TYPE(HashTable_t), POINTER :: Hash
   TYPE(HashEntry_t), POINTER :: Entry

   LOGICAL :: Current

   NULLIFY( Entry )
   IF ( .NOT.ASSOCIATED(Hash) ) RETURN

   Current = .FALSE.
   IF ( ASSOCIATED( Hash % CurrentEntry ) ) THEN
      Current = ASSOCIATED( Hash % CurrentEntry % Next )
   END IF

   IF ( Current ) THEN
      Hash % CurrentEntry => Hash % CurrentEntry % Next
   ELSE
      Hash % CurrentBucket = Hash % CurrentBucket + 1

      DO WHILE(  Hash % CurrentBucket < Hash % BucketSize .AND. &
             .NOT.ASSOCIATED( Hash % Bucket(Hash % CurrentBucket) % Head) )
         Hash % CurrentBucket = Hash % CurrentBucket + 1
      END DO

      IF ( Hash % CurrentBucket >= Hash % BucketSize ) RETURN

      Hash % CurrentEntry => Hash % Bucket(Hash % CurrentBucket) % Head
   END IF

   Entry => Hash % CurrentEntry
 END FUNCTION HashNext

!--------------------------------------------------------------------------
!! Call: void HashStats( HashTable_t *hash )
!!
!! Print info about the hash table organization.
!--------------------------------------------------------------------------
!void HashStats( HashTable_t *hash )
!{
!   HashEntry_t *entry;
!   int *BucketEntries;
!   int n,i,j,MaxEntries=0,MinEntries=1<<30;
!  
!   for( i=0; i<hash->BucketSize; i++ )
!   {
!      n = 0;
!      for( entry = hash->Bucket[i]; entry; entry=entry->next ) n++;
!      MaxEntries = MAX( MaxEntries,n );
!      MinEntries = MIN( MinEntries,n );
!   }
!
!   BucketEntries = (int *)calloc( MaxEntries+1,sizeof(int) );
!
!   for( i=0; i<hash->BucketSize; i++ )
!   {
!      n = 0;
!      for( entry = hash->Bucket[i]; entry; entry=entry->next ) n++;
!      BucketEntries[n]++;
!   }
!
!   fprintf( stdout, "\n\nHash table statistics:\n\n" );
!   fprintf( stdout, "Buckets: % 4d\n", hash->BucketSize );
!   fprintf( stdout, "Entries: % 4d\n", hash->TotalEntries );
!   fprintf( stdout, "Min / Bucket: %d\n", MinEntries );
!   fprintf( stdout, "Max / Bucket: %d\n\n", MaxEntries );
!
!
!   for( n=MinEntries; n<=MaxEntries; n++ )
!      if ( BucketEntries[n] )
!        fprintf( stdout, "% 4d entries in % 4d buckets.\n",
!            n,BucketEntries[n] );
!
!   free( BucketEntries );
!}
END MODULE HashTable


!> \}
